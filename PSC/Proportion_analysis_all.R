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
library(ALDEx2)
library(EnhancedVolcano)
library(VennDiagram)
library(DirichletReg)

# load in dataset
data <- readRDS("/Users/amberbangma/Documents/R/PSC/Data/PSC_processed_oct.rds")
DefaultAssay(data) = "RNA"

# create dimplot with celltypes
DimPlot(data, label = T, repel = T, order = c("T_activated_FOS_low", "Absorptive_enterocyte", "IgA_plasma", "WNT5B+", "Aborptive_TA"), pt.size = 0.5, label.size = 8, raster = F) + NoLegend()
#ggsave("Results/Figures/Dimplot.pdf", width = 33, height = 19)

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
frequencies$state =  apply( frequencies[ ,c(43,44)] , 1 , paste , collapse = "-" )
frequencies = frequencies[,-c(43,44)]

# remove doublets and MT-Hi
frequencies=frequencies[,-c(6,10,15,19,25,32)]

# plot frequencies
df = melt(frequencies)
df$state = factor(df$state, levels = c("HC-NI", "UC-NI", "UC-I", "PSC-NI", "PSC-I"))
df$variable = factor(df$variable, levels = c("Absorptive_enterocyte","Immature_enterocyte","BEST4_enterocyte","DUOX2_enterocyte","Absorptive_TA","PLCG2_TA","Ribo_TA","REG_TA","Cycling_TA","Stem","Enteroendocrine","Immature_goblet","Goblet","Tuft","MAST","glia","WNT2B.","WNT5B.","RSPO3.","Inflammatory_fibroblast","myofibroblasts","pericytes","endothelial","IgG_plasma","IgA_plasma","IgM_plasma","Cycling_B","GC","Follicular_B","CD8T","T_activated_FOS_low","T_activated_FOS_high","CD4_T_memory", "Treg","Cycling_T","APC"))
ggplot(df, aes(x=state, y=value, fill=state)) + 
  geom_violin() +
  geom_boxplot(width=0.1, outlier.size = 0.5) +
  facet_wrap( ~ variable, scales="free") +
  theme_bw(base_size = 14) + 
  scale_fill_brewer(palette="Set1")

#ggsave("Results/Figures/frequencies_violin.pdf", width = 33, height = 19)

# fido's pibble model https://jsilve24.github.io/fido/articles/introduction-to-fido.html
# create counts matrix with state
counts=as.data.frame.matrix(table(data$Final_HTO, data$celltypes))
counts$disease =  sapply(strsplit(rownames(counts),"-"), `[`, 1)
counts$inflammation =  sapply(strsplit(rownames(counts),"-"), `[`, 2)
counts$state =  apply(counts[ ,c(43,44)] , 1 , paste , collapse = "-" )
counts = counts[,-c(43,44)]
counts$state = as.factor(counts$state)
#counts$state = relevel(counts$state, ref = "UC-I") #to set reference, otherwise HC-NI is default

# remove doublets and MT-Hi
counts=counts[,-c(6,10,15,19,25,32)]

# create priors 
X = t(model.matrix(~state, data=counts))
X[,1:5]
Y = t(counts[,-37])
Y[1:5,1:5]
upsilon <- nrow(Y)+3
theta <- matrix(0, nrow(Y)-1, nrow(X))
gamma <- diag(nrow(X))
G <- cbind(diag(nrow(Y)-1), -1)
Xi <- (upsilon-nrow(Y))*G%*%diag(nrow(Y))%*%t(G)
priors <- pibble(Y, X, upsilon, theta, gamma, Xi, n_samples = 3000, step_size = 0.004, max_iter = 20000, b1 = 0.99, eps_f = 1e-13)
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
posterior_summary <- filter(posterior_summary, covariate=="statePSC-NI")
focus <- posterior_summary[sign(posterior_summary$p2.5) == sign(posterior_summary$p97.5),]
focus <- unique(focus$coord)
plot(posterior, par="Lambda", focus.coord = focus, focus.cov = rownames(X)[3]) + geom_vline(xintercept=0, linetype="dashed", color="red")

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

# aldex2
aldex = as_data_frame(t(counts[,-37]))
rownames(aldex) = colnames(counts[,-37])
names = colnames(aldex)
disease = sapply(strsplit(names,"-"), `[`, 1)
inflammation = sapply(strsplit(names,"-"), `[`, 2)
state <- paste(disease, inflammation, sep='-')
x.all <- aldex(aldex[,1:13], state[1:13], mc.samples=16, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)
rownames(x.all) = rownames(aldex)
aldex_PSCI <- filter(x.all, x.all$we.ep <= 0.05)

x.all <- aldex(aldex[,c(1:5,14:29)], state[c(1:5,14:29)], mc.samples=16, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)
rownames(x.all) = rownames(aldex)
aldex_PSCNI <- filter(x.all, x.all$we.ep <= 0.05)

x.all <- aldex(aldex[,c(1:5,30:37)], state[c(1:5,30:37)], mc.samples=16, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)
rownames(x.all) = rownames(aldex)
aldex_UCI <- filter(x.all, x.all$we.ep <= 0.05)

x.all <- aldex(aldex[,c(1:5,38:47)], state[c(1:5,38:47)], mc.samples=16, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)
rownames(x.all) = rownames(aldex)
aldex_UCNI <- filter(x.all, x.all$we.ep <= 0.05)

pibble_PSCI <- state_focus_taxa[["statePSC-I"]]
pibble_PSCI_1 = sapply(strsplit(pibble_PSCI,"_"), `[`, c(2))
pibble_PSCI_2 = sapply(strsplit(pibble_PSCI,"_"), `[`, c(3))
pibble_PSCI <- paste(pibble_PSCI_1, pibble_PSCI_2, sep='_')
rm(pibble_PSCI_2,pibble_PSCI_1)
pibble_PSCNI <- state_focus_taxa[["statePSC-NI"]]
pibble_PSCNI_1 = sapply(strsplit(pibble_PSCNI,"_"), `[`, c(2))
pibble_PSCNI_2 = sapply(strsplit(pibble_PSCNI,"_"), `[`, c(3))
pibble_PSCNI <- paste(pibble_PSCNI_1, pibble_PSCNI_2, sep='_')
rm(pibble_PSCNI_2,pibble_PSCNI_1)
pibble_UCNI <- state_focus_taxa[["stateUC-NI"]]
pibble_UCNI_1 = sapply(strsplit(pibble_UCNI,"_"), `[`, c(2))
pibble_UCNI_2 = sapply(strsplit(pibble_UCNI,"_"), `[`, c(3))
pibble_UCNI <- paste(pibble_UCNI_1, pibble_UCNI_2, sep='_')
rm(pibble_UCNI_2,pibble_UCNI_1)
pibble_UCI <- state_focus_taxa[["stateUC-I"]]
pibble_UCI_1 = sapply(strsplit(pibble_UCI,"_"), `[`, c(2))
pibble_UCI_2 = sapply(strsplit(pibble_UCI,"_"), `[`, c(3))
pibble_UCI <- paste(pibble_UCI_1, pibble_UCI_2, sep='_')
rm(pibble_UCI_2,pibble_UCI_1)

# dirichlet
# create counts matrix with state
counts=as.data.frame.matrix(table(data$Final_HTO, data$celltypes))
counts$disease =  sapply(strsplit(rownames(counts),"-"), `[`, 1)
counts$inflammation =  sapply(strsplit(rownames(counts),"-"), `[`, 2)
counts$state =  apply(counts[ ,c(43,44)] , 1 , paste , collapse = "-" )
counts = counts[,-c(43,44)]
counts$state = as.factor(counts$state)
#counts$state = relevel(counts$state, ref = "UC-I") #to set reference, otherwise HC-NI is default

# remove doublets and MT-Hi
counts=counts[,-c(6,10,15,19,25,32)]

  # Dirichlet multinomial regression to detect changes in cell frequencies
  # formula is not quoted, example: counts ~ condition
  # counts is a [samples x cell types] matrix
  # covariates holds additional data to use in the regression
  #
  # Example:
  # counts = do.call(cbind, tapply(seur@data.info$orig.ident, seur@ident, table))
  # covariates = data.frame(condition=gsub('[12].*', '', rownames(counts)))
  # res = dirichlet_regression(counts, covariates, counts ~ condition)
  
  # Calculate regression
  condition = counts$state
  counts = as.data.frame(counts[,-37])
  Y = DR_data(counts)
  #data = cbind(counts, covariates)
  fit = DirichReg(Y ~ condition)
  
  # Get p-values
  u = summary(fit)
  pvals = u$coef.mat[grep('Intercept', rownames(u$coef.mat), invert=T), 4]
  v = names(pvals)
  pvals = matrix(pvals, ncol=length(u$varnames))
  rownames(pvals) = gsub('condition', '', v[1:nrow(pvals)])
  colnames(pvals) = u$varnames
  fit$pvals = pvals
  fit

  dirichlet <- pvals
  View(dirichlet)