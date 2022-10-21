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
names(rawDat) <- c("length","julian","year","survived","above","diff")
#Get rid of some of the data
x <- rawDat[rawDat$length>=minL & rawDat$length<=maxL,]
x <- x[x$julian>=minJ & x$julian<=maxJ,]



version <- "v11_1"
source("build_cpp_v1.r")

management_dayShift <- -7
management_lengthShift <- 4

n_knots <- 250

version <- "v11.1"
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

myJ <- c(1,0,1,0,0,1) #1 mean estimate, 0 means do not estimate
myL <- c(1,0,1,0,0,1)
myT <- c(1,0,1,1,1,0)
myJL <- c(0,0,0,1,1,0)
myJT <- c(0,0,0,0,0,0)
myJLT <- c(0,0,1,0,1,1)

sdList <- list(est=list(), sd=list())
repList <- list()
optList <- list()
#Simulation flags
#Simulation flags
data$sim <- 0
data$cond_sim <- 0

i <- 3
data$H_flag <- 0
data$j_flag <- myJ[i]
data$l_flag <- myL[i]
data$t_flag <- myT[i]
data$z_jl_flag <- myJL[i]
data$z_jlt_flag <- myJLT[i]
if(myJL[i]>0 | myJLT[i]>0){
  data$H_flag <- 1
}
  
version <- "v11.1"
source(paste0("build_TMB_map_",version,".r"))

random <- c('z_jl', 'z_jlt', 'y_re','j_re', 'l_re')
# 
# version <- "v11_1"
# 
# org_obj <- MakeADFun(data
#                  ,silent = TRUE
#                  , parameters
#                  , random=random
#                  , map = map
#                  , DLL=paste0("spde_aniso_wt_",version)
#                  , sdreport = FALSE)
# org_obj$mesh <- mesh_jl
#   
# opt <- TMBhelper::fit_tmb( org_obj,
#                            loopnum = 1,
#                            newtonsteps = 1,
#                            quiet = TRUE,
#                            getsd = FALSE) # where Obj is a compiled TMB object
# rep <- org_obj$report()

nsim <- 1
simRep_2X2 <- list("noGMRF"=rep(list(rep(list(list()),nsim)),8),
               "GMRF" = rep(list(rep(list(list()),nsim)),8),
               "sim_i" = rep(list(rep(list(list()),nsim)),8))

org_obj$env$data$sim_jre <- rep$j_re
org_obj$env$data$sim_lre <- rep$l_re

sim_data <- org_obj$env$data

#Loop over modeling scenarios
tau2 <- org_obj$env$last.par['log_tau_jl2']
kappa <- org_obj$env$last.par['log_kappa_jl']
H <- org_obj$env$last.par[which(names(org_obj$env$last.par)=='ln_H_input_jl')]

hh <- matrix(c(-0.2,0.2,
               0.2,-0.2,
               2,-2,
               -2,2,
               ),4,2, byrow = TRUE) #aniso
jj <- tau2#log(c(2)) #precision
kk <- kappa#log(c(0.03))


for(i in 1:1){ #nrow(simOut)
  i <- 1
  jcnt <- as.integer(1)
  for(h in 1:nrow(hh)){ #aniso
  for(j in jj){ #precision
    for(k in kk){#range
      set.seed(100*i)
      version <- "v11.1"
      source(paste0("build_TMB_data_pars_",version,".r"))
      
      random <- c('z_jl', 'z_jlt', 'y_re','j_re', 'l_re')
      
        
        #Set the year, day, and length random effects equal to their MLE values for the operating model.
        org_obj$env$last.par['log_tau_jl2'] <- tau2 
        org_obj$env$last.par['log_kappa_jl'] <- kappa
        org_obj$env$last.par[which(names(org_obj$env$last.par)=='ln_H_input_jl')] <- hh[h,]
        sim_i <- org_obj$simulate(complete=TRUE) #simulate data from the original model
        #Run the model without length/day random effects
        sim_data$surv <- sim_i$surv
        
        # sim_data$z_jlt_flag <- 1
        # version <- "v11.1"
        # source(paste0("build_TMB_sim_map_",version,".r"))
        # 
        # version <- "v11_1"
        # sim_obj <- MakeADFun(data=sim_data #create a new TMB model with the simulated data
        #                      , parameters=parameters
        #                      , random=random
        #                      , map = map
        #                      , DLL=paste0("spde_aniso_wt_",version)
        #                      , sdreport = FALSE
        #                      , checkParameterOrder = TRUE
        #                      , silent = TRUE)
        # sim_obj$mesh <- mesh_jl
        # 
        # est <- tryCatch(TMBhelper::fit_tmb(sim_obj, getsd = FALSE, quiet = TRUE),
        #                 error = function(e) e)
        # print(paste("GMRF",jcnt,i,h,j,k))
        # simRep_2X2$GMRF[[(i)]][[(jcnt)]] <- sim_obj$report()



        # sim_data$z_jlt_flag <- 0
        # source(paste0("build_TMB_sim_map_",version,".r"))
        # map <- append(map, list(f_phijt2 = as.factor(NA)
        #                         ,f_phijlt = as.factor(NA)
        #                         ,beta = as.factor(rep(NA,2))
        # )) #No between year correlation
        # sim_obj <- MakeADFun(data=sim_data #create a new TMB model with the simulated data
        #                      , parameters=parameters
        #                      , random=random
        #                      , map = map
        #                      , DLL=paste0("spde_aniso_wt_",version)
        #                      , sdreport = FALSE
        #                      , checkParameterOrder = TRUE
        #                      , silent = TRUE)
        # est <- tryCatch(TMBhelper::fit_tmb(sim_obj, getsd = FALSE, quiet = TRUE),
        #                 error = function(e) e)
        # simRep_2X2$noGMRF[[(i)]][[(jcnt)]] <- sim_obj$report()
        
        simRep_2X2$sim_i[[(i)]][[(jcnt)]] <- sim_i
        jcnt <- (jcnt + 1)
      }
    }
  }
}

source("ggplot_simulated_2D_perturbations.r")
