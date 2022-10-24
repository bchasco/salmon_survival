library(ggplot2)
library(INLA)
library(TMB)
library(DHARMa)
library(mvtnorm)
# library(glmmTMB)
library(tidyr)
library(dplyr)

minL <- 85
maxL <- 125
minJ <- 85
maxJ <- 160

x <- read.csv("C:/NOAA/LARGE_data/chin_survival_1998_2019_subset.csv")

#Rename column to original names
names(x) <- c("length","julian","year","survived","above","diff")
x <- x %>%
  na.omit()
  
#Update the lengths for the bypass fish
x$length[x$above==1] <- x$length[x$above==1] + round(x$diff[x$above==1]/7,0)*2  


#Get rid of some of the data
#Remove fish that are too small or big for the LGR pit tagged fish
x <- x[!(x$above==0 & (x$length<minL | x$length>maxL)),]
#Remove fish that are too small or big for the LGR pit tagged fish
x <- x[!(x$above==1 & (x$length<minL | x$length>=180)),]
#Remove days outside of the min max days
x <- x[!(x$julian<minJ | x$julian>maxJ),]



version <- "v11_5"
source("build_cpp_v1.r")

management_dayShift <- -7
management_lengthShift <- 4

n_knots <- 50

version <- "v11.5"

#1 means use the bypass multivariate
t_bypass <- 1
j_bypass <- 0
l_bypass <- 0
jl_bypass <- 0
AR_flag <- 0

z_flag <- 1
H_flag <- 1
j_flag <- 1
l_flag <- 1
t_flag <- 1
z_jl_flag <- 0
z_jlt_flag <- 0

ad_mar_y <- t_flag
ad_mar_l <- l_flag
ad_mar_j <- j_flag

source(paste0("build_raw_data_",version,".r"))
source(paste0("build_TMB_data_pars_",version,".r"))





#This overrides the bypass argument if flag = FALSE

data$ad_proj <- 1

# myJ <- c(0,0,1,0,0,1) #1 mean estimate, 0 means do not estimate
# myL <- c(0,0,1,0,0,1)
# myT <- c(0,0,1,1,1,0)
# myJL <- c(0,0,0,1,1,0)
# myJT <- c(0,0,0,0,0,0)
# myJLT <- c(1,1,1,0,1,1)
# 
# #Simulation flags
# #Simulation flags
data$sim <- 0
data$cond_sim <- 0
# 
# for(i in c(1:1)){
  # data$H_flag <- 0
  # data$j_flag <- myJ[i]
  # data$l_flag <- myL[i]
  # data$t_flag <- myT[i]
  # data$z_jl_flag <- myJL[i]
  # data$z_jlt_flag <- myJLT[i]
  # if(myJL[i]>0 | myJLT[i]>0){
  #   data$H_flag <- 1
  # }
  
  version <- "v11.5"
  source(paste0("build_TMB_map_",version,".r"))

  random <- c('z_jl', 'z_jlt', 'y_re','j_re', 'l_re')
  
  version <- "v11_5"
  obj <- MakeADFun(data
                   ,silent = FALSE
                   , parameters
                   , random=random
                   , map = map
                   , DLL=paste0("spde_aniso_wt_",version)
                   , sdreport = FALSE)
  obj$mesh <- mesh_jl
  
  # opt <- nlminb(obj$par, obj$fn, obj$gr)
  opt <- TMBhelper::fit_tmb( obj,
                             loopnum = 1,
                             newtonsteps = 1,
                             quiet = TRUE,
                             getsd = FALSE) # where Obj is a compiled TMB object
  rep <- obj$report()
  
  optList[[i]] <- opt
  sd <- sdreport(obj)
  sdList$est[[i]] <- as.list(sd, 'Estimate', report=TRUE)
  sdList$sd[[i]] <- as.list(sd, 'Std. Error', report=TRUE)
  repList[[i]] <- rep
  repList[[i]]$AIC <- 2*opt$objective+2*length(opt$par)
  print(paste("AIC       ", 2*opt$objective+2*length(opt$par)))
  print(paste("gradient       ", obj$gr()))
# }
# save(optList,file = "optList_v11.rData")
# save(sdList,file = "sdList_v11.rData")
# save(repList,file = "repList_v11.rData")
