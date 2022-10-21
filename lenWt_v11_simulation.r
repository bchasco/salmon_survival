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

n_knots <- 250

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

data$ad_mar_y <- data$t_flag
data$ad_mar_l <- data$l_flag
data$ad_mar_j <- data$j_flag
data$ad_proj <- 1

myJ <- c(1,0,1,0,0) #1 mean estimate, 0 means do not estimate
myL <- c(1,0,1,0,0)
myT <- c(1,0,1,1,1)
myJL <- c(0,0,0,1,1)
myJT <- c(0,0,0,0,0)
myJLT <- c(0,0,1,0,1)

# sdList <- list(est=list(), sd=list())
# repList <- list()
# optList <- list()
#Simulation flags
#Simulation flags
data$sim <- 1
data$cond_sim <- 0 #Simulate from all of the distributions

data$H_flag <- 1
data$j_flag <- 1
data$l_flag <- 1
data$t_flag <- 1
data$z_jl_flag <- 0
data$z_jlt_flag <- 1
data$z_jt_flag <- 0
data$H_flag <- 1

source(paste0("build_TMB_map_",version,".r"))
map <- append(map, list(f_phijt2 = as.factor(NA)
                        ,f_phijlt = as.factor(NA)
                        ,beta = as.factor(rep(NA,2))
)) #No between year correlation

random <- c('z_jl', 'z_jlt', 'y_re','j_re', 'l_re', 'z_jt')

#Build the operating model based on the 'true' data 
org_obj <- MakeADFun(data
                 , parameters
                 , random=random 
                 , map = map
                 , DLL=paste0("spde_aniso_wt_",version)
                 , sdreport = FALSE
                 , checkParameterOrder = FALSE)

opt <- TMBhelper::fit_tmb( org_obj,
                           loopnum = 5,
                           newtonsteps = 2,
                           quiet = TRUE,
                           getsd = FALSE) # where Obj is a compiled TMB object
rep <- org_obj$report()

nsim <- 50 
# simRep <- list("noGMRF"=list(),
#                "GMRF" = list(),
#                "sim_i" = list())

sim_data <- data
org_obj$env$data$sim_jre <- rep$j_re
org_obj$env$data$sim_lre <- rep$l_re

#Loop over modeling scenarios
set.seed(10000)
for(i in 1:nsim){ #nrow(simOut)

  source(paste0("build_TMB_data_pars_",version,".r"))
  map <- append(map, list(f_phijt2 = as.factor(NA)
                          ,f_phijlt = as.factor(NA)
                          ,beta = as.factor(rep(NA,2))
  )) #No between year correlation

  
  #Set the year, day, and length random effects equal to their MLE values for the operating model.
  sim_i <- org_obj$simulate(complete=TRUE) #simulate data from the original model
  sim_data$surv <- sim_i$surv
  tmp <- aggregate(list(s=sim_data$surv,t=sim_data$total),by=list(y=sim_data$t_i),sum)
  
  #Run the full model
  sim_data$z_jlt_flag <- 1
  sim_data$z_jl_flag <- 0
  sim_data$z_jt_flag <- 0
  sim_data$t_flag <- 1
  sim_data$l_flag <- 1
  sim_data$j_flag <- 1
  random <- c('z_jl', 'z_jlt', 'y_re','j_re', 'l_re')
  source(paste0("build_TMB_sim_map_",version,".r"))
  map <- append(map, list(f_phijt2 = as.factor(NA)
                          ,f_phijlt = as.factor(NA)
  )) #No between year correlation

  sim_obj <- MakeADFun(data=sim_data #create a new TMB model with the simulated data
                       , parameters=parameters
                       , random=random
                       , map = map
                       , DLL=paste0("spde_aniso_wt_",version)
                       # , sdreport = FALSE
                       , checkParameterOrder = TRUE
                       , silent = TRUE)

  print(paste("Aggregate simulated survival for sim  ", i, round(sum(sim_data$surv)/sum(sim_data$total),4)))
  est <- tryCatch(TMBhelper::fit_tmb(sim_obj, getsd = FALSE, quiet = TRUE),
                  error = function(e) e)
  print(paste(i,2))
  
  if(!is.null(est$objective)){
    rep <- sim_obj$report() #report the simulation
  }
  simRep[[2]][[i]] <- sim_obj$report()

  #Run the model without length/day random effects
  sim_data$z_jlt_flag <- 0
  sim_data$z_jl_flag <- 0
  sim_data$z_jt_flag <- 0
  sim_data$t_flag <- 1
  sim_data$l_flag <- 1
  sim_data$j_flag <- 1
  source(paste0("build_TMB_sim_map_",version,".r"))
  map <- append(map, list(f_phijt2 = as.factor(NA)
                          ,f_phijlt = as.factor(NA)
                          ,beta = as.factor(rep(NA,2))
  )) #No between year correlation
  sim_obj <- MakeADFun(data=sim_data #create a new TMB model with the simulated data
                       , parameters=parameters
                       , random=random
                       , map = map
                       , DLL=paste0("spde_aniso_wt_",version)
                       , sdreport = FALSE
                       , checkParameterOrder = TRUE
                       , silent = TRUE)
  est <- tryCatch(TMBhelper::fit_tmb(sim_obj, getsd = FALSE, quiet = TRUE),
                  error = function(e) e)
  print(paste(i,1))
  
  if(!is.null(est$objective)){
    rep <- sim_obj$report() #report the simulation
  }
  simRep[[1]][[i]] <- sim_obj$report()
  simRep$sim_i[[i]] <- sim_i
}
# 
# save(simRep, file="simRep.rData")
# # par(mfrow=c(3,3))
# # for(i in 1:9)
# #   plot(c(simOut[,,i,2]),c(simOut[,,i,1]))
# # 
# # par(mfrow=c(2,2))
# # plot(c(simDat), c(simOuty[,,1]), main="No GMRF", ylim=c(0,0.3), xlim=c(0,0.3));segments(0,0,1,1)
# # plot(c(simDat), c(simOuty[,,2]), main="GMRF", ylim=c(0,0.3), xlim=c(0,0.3));segments(0,0,1,1)
# # plot(rep(c(simDat),9), c(simOut[,,,1]), main="No GMRF", ylim=c(0,0.3), xlim=c(0,0.3));segments(0,0,1,1)
# # plot(rep(c(simDat),9), c(simOut[,,,2]), main="GMRF", ylim=c(0,0.3), xlim=c(0,0.3));segments(0,0,1,1)
# # 
