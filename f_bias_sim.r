f_bias_sim <- function(obj = obj,
                           TMB_list = TMB_list,
                           map = map,
                           bias_sim_n = bias_sim_n){
  
  sim_out <- matrix(NA,bias_sim_n,length(obj$par))

  for(i in 1:bias_sim_n){
    #Set the year, day, and length random effects equal to their MLE values for the operating model.
    sim_i <- obj$simulate(complete=TRUE) #simulate data from the original model
    print(paste(sum(sim_i$surv), "sim ", i , "out of ", bias_sim_n))
    random <- obj$env$.random
    sim_obj <- MakeADFun(data=sim_i #create a new TMB model with the simulated data
                         , parameters=TMB_list$parameters
                         , random=random
                         , map = map
                         , DLL=paste0("spde_aniso_wt_","v11_6")
                         # , sdreport = FALSE
                         , checkParameterOrder = TRUE
                         , silent = TRUE)
    opt <- TMBhelper::fit_tmb( sim_obj,
                               loopnum = 1,
                               newtonsteps = 1,
                               quiet = TRUE,
                               getsd = FALSE) # where Obj is a compiled TMB object
    
    sim_out[i,] <- opt$par
  }
  return(sim_out)
}

