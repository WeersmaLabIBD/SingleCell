#subset IEL, LPL and BLOOD cells

#subset allpts_allcells for 1. colname contains BLOOD, 2. colname contains IEL, 3. colname contains LPL

library(dplyr)
setwd("..")
load("fourth_file.Robj")

library(dplyr)
IEL_allpts<-select(fourth_file, contains("IEL"))

LPL_allpts<-select(fourth_file, contains("LPL"))
  
BLOOD_allpts<-select(fourth_file, contains("BLOOD"))  

