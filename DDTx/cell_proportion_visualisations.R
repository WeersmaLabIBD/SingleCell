######
cellnumbers<-read.table("~/xx/cell_numbers_compartment.tsv", header=T)
cellnumbers$proportion_recipient_immune_cells2=as.numeric(as.character(gsub("%","",cellnumbers$proportion_recipient_immune_cells)
                                                                       + ))

cellnumbers$my_order[cellnumbers$ACR_grade=="Zero"]=0
cellnumbers$my_order[cellnumbers$ACR_grade=="Zero_one"]=1
cellnumbers$my_order[cellnumbers$ACR_grade=="Zero_one"]=0.5
cellnumbers$my_order[cellnumbers$ACR_grade=="One"]=1
cellnumbers$my_order[cellnumbers$ACR_grade=="Two"]=2
cellnumbers$my_order[cellnumbers$ACR_grade=="Three"]=3
cellnumbers$my_order[cellnumbers$ACR_grade=="Post-three"]=3.5

# ACR category immune cell proportions recipient cells of all immune cells
######
ggplot(cellnumbers, aes(x = as.factor(my_order),y = proportion_recipient_immune_cells2)) +
  geom_boxplot() + theme_bw() + xlab ("ACR grade") + ylab ("% recipient immune  cells") + geom_smooth(method='lm',formula= y~x)
#####
ggplot(cellnumbers, aes(y = proportion_recipient_immune_cells2, x = timepoint, color = as.factor(my_order))) + geom_point() + theme_bw()

ggplot(cellnumbers, aes(y = proportion_recipient_immune_cells2, x = timepoint, color = patient)) + geom_point() + theme_bw()


cellnumbers_celltypes<-read.table("~/xx/cell_numbers_adult_elmentaite_martin_immune_analysis.txt", header=T, dec=",")


ggplot(cellnumbers_celltypes, aes(y =Proportion_immune_of_all, x = timepoint, color = donor_recipient)) + geom_point() + theme_bw()


celltypes<-c("Activated_CD4_T", "Activated_CD8_T", "Adult_Glia","BEST4_positive_epithelial", "cDC1", "cDC2","IgA_plasma_cell","IgG_plasma_cell","ILC3","LYVE1_positive_Macrophage","Macrophages","Memory_B","Monocytes", "Naive_B","TA" ,"TRGV2_gdT","Enterocyte","arterial_capillary","Cycling_B_cell","Mast_cell","Paneth","Contractile_pericyte_PLN_positive","D_cells_SST_positive","Goblet_cell","I_cells_CCK_positive","K_cells_GIP_positive","LEC1_ACKR4_positive","Stem_cells","Tuft","Mature_arterial_EC" ,"Mature_venous_EC" ,"myofibroblast_RSPO2_positive","Stromal_1_ADAMDEC1_positive","Stromal_1_CCL11_positive","Stromal_2_NPY_positive","Stromal_3_C7_positive","Transitional_Stromal_3_C3_positive","LEC3_ADGRG3_positive","gdT","TRGV5_7_gdT","CX3CR1_positive_CD8_Tmem","TRGV4_gdT")  
t_celltypes<-c("Activated_CD4_T", "Activated_CD8_T","TRGV2_gdT","gdT","TRGV5_7_gdT","CX3CR1_positive_CD8_Tmem","TRGV4_gdT")  
epithelial_celltypes<-c("Adult_Glia","BEST4_positive_epithelial","TA" ,"Enterocyte","Paneth","D_cells_SST_positive","I_cells_CCK_positive","K_cells_GIP_positive","Goblet_cell","Stem_cells","Tuft")
myeloid_celltypes<-c("cDC1", "cDC2","ILC3","LYVE1_positive_Macrophage","Macrophages","Monocytes","Mast_cell")
b_celltypes<-c("IgA_plasma_cell","IgG_plasma_cell","Memory_B","Naive_B","Cycling_B_cell")
stromal_celltypes<-c("arterial_capillary","Contractile_pericyte_PLN_positive","LEC1_ACKR4_positive","Mature_arterial_EC" ,"Mature_venous_EC" ,"myofibroblast_RSPO2_positive","Stromal_1_ADAMDEC1_positive","Stromal_1_CCL11_positive","Stromal_2_NPY_positive","Stromal_3_C7_positive","Transitional_Stromal_3_C3_positive","LEC3_ADGRG3_positive")
immune_celltypes<-c("Activated_CD4_T", "Activated_CD8_T","TRGV2_gdT","gdT","TRGV5_7_gdT","CX3CR1_positive_CD8_Tmem","TRGV4_gdT", "cDC1", "cDC2","ILC3","LYVE1_positive_Macrophage","Macrophages","Monocytes","Mast_cell","IgA_plasma_cell","IgG_plasma_cell","Memory_B","Naive_B","Cycling_B_cell" )

colnames(cellnumbers_celltypes)
prop_categories<-cellnumbers_celltypes[c(31,41,58,59,54:56)] 
y<-prop_categories
library(reshape2)
y<-melt(y, id.vars=c("donor_recipient","ACR_grade", "sample", "timepoint" ))
y[y == "NaN"] <- "0"
y[y == "Inf"] <- "0"

y$value<-as.numeric(y$value)
y_pt1<-y[y$sample == "Pt1",]
y_pt2<-y[y$sample == "Pt2",]
y_pt3<-y[y$sample == "Pt3",]
y_pt1_d<-y_pt1[y_pt1$donor_recipient == "d",]
y_pt1_d$value<-as.numeric(y_pt1_d$value)

y_pt2_d<-y_pt2[y_pt2$donor_recipient == "d",]
y_pt2_d$value<-as.numeric(y_pt2_d$value)

y_pt3_d<-y_pt3[y_pt3$donor_recipient == "d",]
y_pt3_d$value<-as.numeric(y_pt3_d$value)


y_pt1_r<-y_pt1[y_pt1$donor_recipient == "r",]
y_pt1_r$value<-as.numeric(y_pt1_r$value)

y_pt2_r<-y_pt2[y_pt2$donor_recipient == "r",]
y_pt2_r$value<-as.numeric(y_pt2_r$value)

y_pt3_r<-y_pt3[y_pt3$donor_recipient == "r",]
y_pt3_r$value<-as.numeric(y_pt3_r$value)

library(ggplot2)
ggplot(y_pt1_d, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 

ggplot(y_pt2_d, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 

ggplot(y_pt3_d, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 


ggplot(y_pt1_r, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 

ggplot(y_pt2_r, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 

ggplot(y_pt3_r, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 

ggplot(y_pt3_r, aes(x=ACR_grade,y=value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 


prop_celltypes<-cellnumbers_celltypes

for(i in immune_celltypes){
  prop_celltypes[, ncol(prop_celltypes) + 1] <- prop_celltypes[i]/prop_celltypes$Total_immune
  names(prop_celltypes)[ncol(prop_celltypes)] <- paste0("prop_", i)
}
prop_celltypes
colnames(prop_celltypes)

prop_celltypes<-prop_celltypes[c(31,41,58,59,60:78)]

x<-prop_celltypes

library(reshape2)
x<-melt(x, id.vars=c("donor_recipient","ACR_grade", "sample", "timepoint" ))
x[x == "NaN"] <- "0"
x[x == "Inf"] <- "0"


library(dplyr)
x$value<-as.numeric(x$value)
x_pt1<-x[x$sample == "Pt1",]
x_pt2<-x[x$sample == "Pt2",]
x_pt3<-x[x$sample == "Pt3",]
x_pt1_d<-x_pt1[x_pt1$donor_recipient == "d",]
x_pt1_d$value<-as.numeric(x_pt1_d$value)

x_pt2_d<-x_pt2[x_pt2$donor_recipient == "d",]
x_pt2_d$value<-as.numeric(x_pt2_d$value)

x_pt3_d<-x_pt3[x_pt3$donor_recipient == "d",]
x_pt3_d$value<-as.numeric(x_pt3_d$value)


x_pt1_r<-x_pt1[x_pt1$donor_recipient == "r",]
x_pt1_r$value<-as.numeric(x_pt1_r$value)

x_pt2_r<-x_pt2[x_pt2$donor_recipient == "r",]
x_pt2_r$value<-as.numeric(x_pt2_r$value)

x_pt3_r<-x_pt3[x_pt3$donor_recipient == "r",]
x_pt3_r$value<-as.numeric(x_pt3_r$value)

#celltype  proportions
library(ggplot2)
ggplot(x_pt1_d, aes(x=timepoint,value, fill=variable))+
geom_bar(position=position_stack(), stat="identity") 
ggplot(x_pt1_d, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(x_pt2_d, aes(x=timepoint,value, fill=variable))+
geom_bar(position=position_stack(), stat="identity") 
ggplot(x_pt2_d, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(x_pt3_d, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 
ggplot(x_pt3_d, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(x_pt1_r, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 
z<-x_pt1_r[x_pt1_r$timepoint != c(47,121),]
ggplot(z, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(x_pt2_r, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 
z<-x_pt2_r[x_pt2_r$timepoint != 12,]
ggplot(z, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(x_pt3_r, aes(x=timepoint,value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity") 

#celltype category proportions
#####
ggplot(y_pt1_d, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(y_pt2_d, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(y_pt3_d, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()


ggplot(y_pt1_r, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()
z<-y_pt1_r[y_pt1_r$timepoint != c(47,121),]
ggplot(z, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(y_pt2_r, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()
z<-y_pt2_r[y_pt2_r$timepoint != 12,]
ggplot(z, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()


ggplot(y_pt3_r, aes(x=timepoint,value, fill=variable))+
  geom_area(position=position_stack(), stat="identity", color="black") + theme_bw()

ggplot(y_pt3_r, aes(x=ACR_grade,y=value, fill=variable))+
  geom_bar(position=position_stack(), stat="identity", color="black") + theme_bw()

#cell types per ACR grade
######

x$my_order[x$ACR_grade=="Zero"]=0
x$my_order[x$ACR_grade=="Zero_one"]=1
x$my_order[x$ACR_grade=="Zero_one"]=0.5
x$my_order[x$ACR_grade=="One"]=1
x$my_order[x$ACR_grade=="Two"]=2
x$my_order[x$ACR_grade=="Three"]=3
x$my_order[x$ACR_grade=="Post-three"]=3.5

x_r<-x[x$donor_recipient == "r",]
x_d<-x[x$donor_recipient == "d",]

x_d$ACR_prop<-x_d$value
x_d$ACR_prop[x_d$ACR_grade=="Zero"]<-x_d$value[x_d$ACR_grade=="Zero"]/8
x_d$ACR_prop[x_d$ACR_grade=="Zero_one"]=x_d$value[x_d$ACR_grade=="Zero_one"]/5
x_d$ACR_prop[x_d$ACR_grade=="One"]=x_d$value[x_d$ACR_grade=="One"]/8
x_d$ACR_prop[x_d$ACR_grade=="Two"]=x_d$value[x_d$ACR_grade=="Two"]/2
x_d$ACR_prop[x_d$ACR_grade=="Three"]=x_d$value[x_d$ACR_grade=="Three"]/2
x_d$ACR_prop[x_d$ACR_grade=="Post-three"]=x_d$value[x_d$ACR_grade=="Post-three"]/3

ggplot(x_d, aes(x=my_order,y=ACR_prop, fill=variable))+
  geom_bar(position=position_stack(), stat="identity", color="black") + theme_bw()


x_r$ACR_prop<-x_r$value
x_r$ACR_prop[x_r$ACR_grade=="Zero"]<-x_r$value[x_r$ACR_grade=="Zero"]/2
x_r$ACR_prop[x_r$ACR_grade=="Zero_one"]=x_r$value[x_r$ACR_grade=="Zero_one"]/5
x_r$ACR_prop[x_r$ACR_grade=="One"]=x_r$value[x_r$ACR_grade=="One"]/3
x_r$ACR_prop[x_r$ACR_grade=="Two"]=x_r$value[x_r$ACR_grade=="Two"]/2
x_r$ACR_prop[x_r$ACR_grade=="Three"]=x_r$value[x_r$ACR_grade=="Three"]/2
x_r$ACR_prop[x_r$ACR_grade=="Post-three"]=x_r$value[x_r$ACR_grade=="Post-three"]/4

ggplot(x_r, aes(x=my_order,y=ACR_prop, fill=variable))+
  geom_bar(position=position_stack(), stat="identity", color="black") + theme_bw()


ggplot(x_d, aes(x=my_order,y=value, fill=variable))
  geom_bar(position = position_stack(), stat="identity", color="black") + theme_bw()

##
  y$my_order[y$ACR_grade=="Zero"]=0
  y$my_order[y$ACR_grade=="Zero_one"]=1
  y$my_order[y$ACR_grade=="Zero_one"]=0.5
  y$my_order[y$ACR_grade=="One"]=1
  y$my_order[y$ACR_grade=="Two"]=2
  y$my_order[y$ACR_grade=="Three"]=3
  y$my_order[y$ACR_grade=="Post-three"]=3.5
  
  y_r<-y[y$donor_recipient == "r",]
  y_d<-y[y$donor_recipient == "d",]
  
  y_d$ACR_prop<-y_d$value
  y_d$ACR_prop[y_d$ACR_grade=="Zero"]<-y_d$value[y_d$ACR_grade=="Zero"]/8
  y_d$ACR_prop[y_d$ACR_grade=="Zero_one"]=y_d$value[y_d$ACR_grade=="Zero_one"]/5
  y_d$ACR_prop[y_d$ACR_grade=="One"]=y_d$value[y_d$ACR_grade=="One"]/8
  y_d$ACR_prop[y_d$ACR_grade=="Two"]=y_d$value[y_d$ACR_grade=="Two"]/2
  y_d$ACR_prop[y_d$ACR_grade=="Three"]=y_d$value[y_d$ACR_grade=="Three"]/2
  y_d$ACR_prop[y_d$ACR_grade=="Post-three"]=y_d$value[y_d$ACR_grade=="Post-three"]/3
  
  ggplot(y_d, aes(x=my_order,y=ACR_prop, fill=variable))+
    geom_bar(position=position_stack(), stat="identity", color="black") + theme_bw()
  
  
  y_r$ACR_prop<-y_r$value
  y_r$ACR_prop[y_r$ACR_grade=="Zero"]<-y_r$value[y_r$ACR_grade=="Zero"]/2
  y_r$ACR_prop[y_r$ACR_grade=="Zero_one"]=y_r$value[y_r$ACR_grade=="Zero_one"]/5
  y_r$ACR_prop[y_r$ACR_grade=="One"]=y_r$value[y_r$ACR_grade=="One"]/3
  y_r$ACR_prop[y_r$ACR_grade=="Two"]=y_r$value[y_r$ACR_grade=="Two"]/2
  y_r$ACR_prop[y_r$ACR_grade=="Three"]=y_r$value[y_r$ACR_grade=="Three"]/2
  y_r$ACR_prop[y_r$ACR_grade=="Post-three"]=y_r$value[y_r$ACR_grade=="Post-three"]/4
  
  ggplot(y_r, aes(x=my_order,y=ACR_prop, fill=variable))+
    geom_bar(position=position_stack(), stat="identity", color="black") + theme_bw()
  
a<-y_r  
  
a$pt_tp_dr<-paste(y_r$sample, y_r$timepoint, y_r$donor_recipient, sep="_")

ggplot(a, aes(x=my_order,y=ACR_prop, fill=variable))+
  geom_bar(position = position_dodge(), stat="identity", color="black") + theme_bw()


