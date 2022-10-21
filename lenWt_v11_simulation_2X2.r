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

# #Build the operating model based on the 'true' data
org_obj <- MakeADFun(data
                 , parameters
                 , random=random
                 , map = map
                 , DLL=paste0("spde_aniso_wt_",version)
                 , sdreport = FALSE
                 , checkParameterOrder = FALSE)

opt <- TMBhelper::fit_tmb( org_obj,
                           loopnum = 1,
                           newtonsteps = 1,
                           quiet = TRUE,
                           getsd = FALSE) # where Obj is a compiled TMB object
rep <- org_obj$report()

nsim <- 2
simRep_2X2 <- list("noGMRF"=rep(list(rep(list(list()),nsim)),8),
               "GMRF" = rep(list(rep(list(list()),nsim)),8),
               "sim_i" = rep(list(rep(list(list()),nsim)),8))

sim_data <- data
org_obj$env$data$sim_jre <- rep$j_re
org_obj$env$data$sim_lre <- rep$l_re

#Loop over modeling scenarios
tau2 <- org_obj$env$last.par['log_tau2_jl']
kappa <- org_obj$env$last.par['log_kappa_jl']
H <- org_obj$env$last.par[which(names(org_obj$env$last.par)=='ln_H_input_jl')]

hh <- c(-0.5,-1) #aniso
jj <- log(c(2)) #precision
kk <- log(c(0.03))

for(i in 1:1){ #nrow(simOut)
  jcnt <- as.integer(1)
  for(h in 1:1){ #aniso
  for(j in jj){ #precision
    for(k in kk){#range
      set.seed(100*i)
      source(paste0("build_TMB_data_pars_",version,".r"))
        map <- append(map, list(f_phijt2 = as.factor(NA)
                                ,f_phijlt = as.factor(NA)
                                ,beta = as.factor(rep(NA,2))
        )) #No between year correlation
        
        
        #Set the year, day, and length random effects equal to their MLE values for the operating model.
        org_obj$env$last.par['log_tau2_jl'] <- tau2 
        org_obj$env$last.par['log_kappa_jl'] <- kappa
        org_obj$env$last.par[which(names(org_obj$env$last.par)=='ln_H_input_jl')] <- hh
        sim_i <- org_obj$simulate(complete=TRUE) #simulate data from the original model
        #Run the model without length/day random effects
        sim_data$surv <- sim_i$surv
        
        sim_data$z_jlt_flag <- 1
        source(paste0("build_TMB_sim_map_",version,".r"))
        map <- append(map, list(f_phijt2 = as.factor(NA)
                                ,f_phijlt = as.factor(NA)
                                ,beta = as.factor(rep(NA,2))
        )) #No between year correlation
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
        print(paste("GMRF",jcnt,i,h,j,k))
        simRep_2X2$GMRF[[(i)]][[(jcnt)]] <- sim_obj$report()
        

        
        sim_data$z_jlt_flag <- 0
        source(paste0("build_TMB_sim_map_",version,".r"))
        map <- append(map, list(f_phijt2 = as.factor(NA)
                                ,f_phijlt = as.factor(NA)
                                ,beta = as.factor(rep(NA,2))
        )) #No between year correlation
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

        simRep_2X2$noGMRF[[(i)]][[(jcnt)]] <- sim_obj$report()
        simRep_2X2$sim_i[[(i)]][[(jcnt)]] <- sim_i
        jcnt <- (jcnt + 1)
      }
    }
  }
}
