library(ggplot2)
library(INLA)
library(TMB)
library(DHARMa)
library(mvtnorm)

minL <- 85
maxL <- 125
minJ <- 85
maxJ <- 160

if(!exists("rawDat")){ #This is a huge file only read it if it's necessary
  rawDat <- read.csv("master_tag_with_weight.csv")
}

#Get rid of some of the data
x <- rawDat[rawDat$length>=minL & rawDat$length<=maxL,]
x <- x[x$julian>=minJ & x$julian<=maxJ,]
x <- x[x$trans=="river",]

version <- "v11"
source("build_cpp_v1.r")

management_dayShift <- -7
management_lengthShift <- 4

source(paste0("build_raw_data_",version,".r"))
source(paste0("build_TMB_data_pars_",version,".r"))


data$AR_flag <- 1
data$z_flag <- 1
data$H_flag <- 0

data$j_flag <- 1
data$l_flag <- 1
data$t_flag <- 1
data$z_l_flag <- 1
data$z_j_flag <- 1
data$z_jl_flag <- 0
data$z_jlt_flag <- 0
data$z_jt_flag <- 0

data$ad_mar_y <- 0
data$ad_mar_l <- 0
data$ad_mar_j <- 0
data$ad_proj <- 0

myJ <- c(1,1,1,0,0) #1 mean estimate, 0 means do not estimate
myL <- c(1,1,1,0,0)
myT <- c(1,1,1,1,1)
myJL <- c(0,0,0,1,1)
myJT <- c(0,1,0,0,0)
myJLT <- c(0,0,1,0,1)

sdList <- list(est=list(), sd=list())
repList <- list()
optList <- list()
#Simulation flags
data$sim <- 1
data$cond_sim <- 1

for(i in 3:3){
  data$H_flag <- 0
  data$j_flag <- myJ[i]
  data$l_flag <- myL[i]
  data$t_flag <- myT[i]
  data$z_jl_flag <- myJL[i]
  data$z_jlt_flag <- myJLT[i]
  data$z_jt_flag <- myJT[i]
  if(myJL[i]>0 | myJLT[i]>0){
    data$H_flag <- 1
  }
  source(paste0("build_TMB_map_",version,".r"))
  map <- append(map, list(f_phijt2 = as.factor(NA)
                            ,f_phijlt = as.factor(NA)
                          )) #No between year correlation
  
  random <- c('z_jl', 'z_jlt', 'y_re','j_re', 'l_re', 'z_jt')
  
  obj <- MakeADFun(data
                   # ,silent = TRUE
                   , parameters
                   , random=random 
                   , map = map
                   , DLL=paste0("spde_aniso_wt_",version)
                   , sdreport = FALSE)
  
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  rep <- obj$report()
  
  sim <- replicate(1, {
    simdata <- obj$simulate()
    simdata <- simdata$surv
  })

  d_res <- createDHARMa(sim[,],
                        data$surv,
                        fittedPredictedResponse = plogis(rep$eta_i)*data$total)

  png("lenWt_v11_Dharama_sim.png", unit="in", height=5, width = 5, res = 600)
  KSi <- testUniformity(d_res)
  dev.off()
}
