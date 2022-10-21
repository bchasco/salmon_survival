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

# if(!exists("rawDat")){ #This is a huge file only read it if it's necessary
  rawDat <- read.csv("C:/NOAA/LARGE_data/chin_survival_1998_2019_subset.csv")
  # rawDat <- read.csv("C:/NOAA/LARGE_data/doublecount.csv")
  # }
# rawDat <- rawDat[!(rawDat$Migration.Year.YYYY!=2001 & rawDat$tagged_above_LGR==1),] 

#Rename column to original names
names(rawDat) <- c("length","julian","year","survived","above","diff")
#Get rid of some of the data
x <- rawDat[rawDat$length>=minL & rawDat$length<=maxL,]
x <- x[x$julian>=minJ & x$julian<=maxJ,]



version <- "v11_5"
source("build_cpp_v1.r")

management_dayShift <- -7
management_lengthShift <- 4

n_knots <- 250

version <- "v11.5"
source(paste0("build_raw_data_",version,".r"))
source(paste0("build_TMB_data_pars_",version,".r"))


data$AR_flag <- 1
data$z_flag <- 1
data$H_flag <- 1

data$j_flag <- 1
data$l_flag <- 1
data$t_flag <- 1
data$z_l_flag <- 1
data$z_j_flag <- 1
data$z_jl_flag <- 0
data$z_jlt_flag <- 0
data$z_jt_flag <- 0

data$ad_mar_y <- data$t_flag
data$ad_mar_l <- data$l_flag
data$ad_mar_j <- data$j_flag
data$ad_proj <- 1

myJ <- c(1,0,1,0,0,1) #1 mean estimate, 0 means do not estimate
myL <- c(1,0,1,0,0,1)
myT <- c(1,0,1,1,1,0)
myJL <- c(0,0,0,1,1,0)
myJT <- c(0,0,0,0,0,0)
myJLT <- c(0,1,1,0,1,1)

sdList <- list(est=list(), sd=list())
repList <- list()
optList <- list()
#Simulation flags
#Simulation flags
data$sim <- 0
data$cond_sim <- 0

for(i in c(3:3)){
  data$H_flag <- 0
  data$j_flag <- myJ[i]
  data$l_flag <- myL[i]
  data$t_flag <- myT[i]
  data$z_jl_flag <- myJL[i]
  data$z_jlt_flag <- myJLT[i]
  if(myJL[i]>0 | myJLT[i]>0){
    data$H_flag <- 1
  }
  
  version <- "v11.3"
  source(paste0("build_TMB_map_",version,".r"))

  random <- c('z_jl', 'z_jlt', 'y_re','j_re', 'l_re')
  
  version <- "v11_5"
  obj <- MakeADFun(data
                   ,silent = TRUE
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
}
# save(optList,file = "optList_v11.rData")
# save(sdList,file = "sdList_v11.rData")
# save(repList,file = "repList_v11.rData")
