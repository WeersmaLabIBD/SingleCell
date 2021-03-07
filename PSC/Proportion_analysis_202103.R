library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)
library(reshape2)
library(phyloseq)
library(fido)
library(compositions)
library(corrplot)

# load in dataset
data<- readRDS("PSC_processed_march_2021.rds")
DefaultAssay(data) = "RNA"

# extract counts to matrix
counts=as.data.frame.matrix(table(data$Final_HTO, data$celltypes))

# calculate proportions
total_row = apply(counts, 1, sum)
pcts = lapply(counts, function(x) {
  x / total_row
})
frequencies = as.data.frame(pcts)
rownames(counts) = rownames(frequencies)

# add state
frequencies$disease =  sapply(strsplit(rownames(frequencies),"-"), `[`, 1)
frequencies$inflammation =  sapply(strsplit(rownames(frequencies),"-"), `[`, 2)
frequencies$state =  apply( frequencies[ ,c(44,45)] , 1 , paste , collapse = "-" )
frequencies = frequencies[,-c(44,45)]

# remove doublets and MT-Hi
frequencies=frequencies[,-c(6,14,20,26,33)]

# plot frequencies
list_celltypes<-list(colnames(frequencies[-39]))
df = melt(frequencies)
df$state = factor(df$state, levels = c("HC-NI", "UC-NI", "UC-I", "PSC-NI", "PSC-I"))
df$variable = factor(df$variable, levels = list_celltypes[[1]])
ggplot(df, aes(x=state, y=value, fill=state)) + 
  geom_boxplot(outlier.size = 1) + 
  facet_wrap( ~ variable, scales="free") +
  theme_bw(base_size = 14) + 
  scale_fill_brewer(palette="Set1")

ggsave("Results/Figures/frequencies.pdf", width = 33, height = 19)

# fido's pibble model https://jsilve24.github.io/fido/articles/introduction-to-fido.html
# create counts matrix with state
counts=as.data.frame.matrix(table(data$Final_HTO, data$celltypes))
counts$disease =  sapply(strsplit(rownames(counts),"-"), `[`, 1)
counts$inflammation =  sapply(strsplit(rownames(counts),"-"), `[`, 2)
counts$state =  apply(counts[ ,c(44,45)] , 1 , paste , collapse = "-" )
counts = counts[,-c(44,45)]
counts$state = as.factor(counts$state)
#counts$state = relevel(counts$state, ref = "UC-I") #to set reference, otherwise HC-NI is default

# remove doublets and MT-Hi
counts=counts[,-c(6,14,20,26,33)]

# create priors 
X = t(model.matrix(~state, data=counts))
X[,1:5]
Y = t(counts[,-39])
Y[1:5,1:5]
upsilon <- nrow(Y)+3
theta <- matrix(0, nrow(Y)-1, nrow(X))
gamma <- diag(nrow(X))
G <- cbind(diag(nrow(Y)-1), -1)
Xi <- (upsilon-nrow(Y))*G%*%diag(nrow(Y))%*%t(G)
priors <- pibble(Y, X, upsilon, theta, gamma, Xi)
print(priors)

# change model to clr
priors <- to_clr(priors)
summary(priors, pars="Lambda")
names_covariates(priors) <- rownames(X)

#plot log probabilites for each state
plot(priors, par="Lambda") + ggplot2::xlim(c(-10, 10))

#check if model fits data (checked by Johannes)
priors$Y <- Y
posterior <- refit(priors, optim_method="adam")
ppc(posterior) + ggplot2::coord_cartesian(ylim=c(0, 500)) #I changed axis from 30000 to 500 for visualization
ppc_summary(posterior)
ppc(posterior, from_scratch=TRUE) +ggplot2::coord_cartesian(ylim=c(0, 500))
ppc_summary(posterior, from_scratch=TRUE)

# look at posterior distribution of regression parameters
posterior_summary <- summary(posterior, pars="Lambda")$Lambda

# focus on coordinates with non-zero effect ("significant") for UC-I 
posterior_summary <- filter(posterior_summary, covariate=="stateUC-I")
focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
focus <- unique(focus$coord)
plot(posterior, par="Lambda", focus.coord = focus, focus.cov = rownames(X)[4]) + geom_vline(xintercept=0, linetype="dashed", color="red")

# create list of significant celltypes for each state
posterior_summary <- summary(posterior, pars="Lambda")$Lambda
state_focus_taxa <- vector("list", length(rownames(X)[2:5]))
names(state_focus_taxa) <- rownames(X)[2:5]
for(state in rownames(X)[2:5]) {
  post_summary <- filter(posterior_summary, covariate %in% state)
  focus <- post_summary[sign(post_summary$p2.5) == sign(post_summary$p97.5),]
  focus <- unique(focus$coord)
  state_focus_taxa[[state]] <- focus
}

# create square plot for differences
mat_state <- diag(length(rownames(X)[2:5]))
rownames(mat_state) <- colnames(mat_state) <- rownames(X)[2:5]
state_pairs <- data.frame(t(combn(rownames(X)[2:5], 2)))
for(pair in 1:nrow(state_pairs)){
  print(pair)
  mat_state[state_pairs[pair,]$X1, state_pairs[pair,]$X2] <- length(setdiff(state_focus_taxa[[state_pairs[pair,]$X1]], state_focus_taxa[[state_pairs[pair,]$X2]]))
  mat_state[state_pairs[pair,]$X2, state_pairs[pair,]$X1] <- length(setdiff(state_focus_taxa[[state_pairs[pair,]$X2]], state_focus_taxa[[state_pairs[pair,]$X1]]))
}
corrplot::corrplot(mat_state, method="square", is.corr=F, addCoef.col = "black", diag = F, tl.col="black", cl.pos="n")

# create square plot for agreement
mat_state <- diag(length(rownames(X)[2:5]))
rownames(mat_state) <- colnames(mat_state) <- rownames(X)[2:5]
state_pairs <- data.frame(t(combn(rownames(X)[2:5], 2)))
for(pair in 1:nrow(state_pairs)){
  print(pair)
  mat_state[state_pairs[pair,]$X1, state_pairs[pair,]$X2] <- length(intersect(state_focus_taxa[[state_pairs[pair,]$X1]], state_focus_taxa[[state_pairs[pair,]$X2]]))
  mat_state[state_pairs[pair,]$X2, state_pairs[pair,]$X1] <- length(intersect(state_focus_taxa[[state_pairs[pair,]$X2]], state_focus_taxa[[state_pairs[pair,]$X1]]))
}
corrplot::corrplot(mat_state, method="square", type = "lower", is.corr=F, addCoef.col = "black", diag = F, tl.col="black", cl.pos="n")

