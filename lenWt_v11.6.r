library(ggplot2)
library(INLA)
library(TMB)
library(DHARMa)
library(mvtnorm)
# library(glmmTMB)

minL <- 85
maxL <- 125
minJ <- 85
maxJ <- 160

rawDat <- read.csv("C:/NOAA/LARGE_data/chin_survival_1998_2019_subset.csv")

#Rename column to original names
names(rawDat) <- c("length","julian","year","survived","above","diff")
#Get rid of some of the data
x <- rawDat[rawDat$length>=minL & rawDat$length<=maxL,]
x <- x[x$julian>=minJ & x$julian<=maxJ,]



version <- "v11_6"
source("build_cpp_v1.r")

n_knots <- 250
version <- "v11.6"
source(paste0("build_raw_data_",version,".r"))
source(paste0("build_TMB_data_pars_",version,".r"))



data$j_flag <- 0
data$l_flag <- 0
data$t_flag <- 0
data$j_bypass <- 1
data$l_bypass <- 0
data$t_bypass <- 0
source(paste0("build_TMB_map_",version,".r"))

random <- c('y_re','j_re', 'l_re')
  
version <- "v11_6"
obj <- MakeADFun(data
                 ,silent = TRUE
                 , parameters
                 , random=random
                 , map = map
                 , DLL=paste0("spde_aniso_wt_",version)
                 , sdreport = FALSE)
# obj$mesh <- mesh_jl
opt <- TMBhelper::fit_tmb( obj,
                           loopnum = 1,
                           newtonsteps = 1,
                           # quiet = TRUE,
                           getsd = FALSE) # where Obj is a compiled TMB object
 rep <- obj$report()
 