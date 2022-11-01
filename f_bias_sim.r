f_bias_sim <- function(obj = obj,
                           TMB_list = TMB_list,
                           map = map,
                           bias_sim_n = bias_sim_n,
                       sim_size = sim_size){
  
  
  sim_out <- matrix(NA,bias_sim_n,length(obj$par))
  sim_proj <- array(NA,c(bias_sim_n,dim(obj$rep$proj_y)))
  
  obj$env$data$sim_size <- sim_size
  
  for(i in 1:bias_sim_n){
    sys <- Sys.time()
    #Set the year, day, and length random effects equal to their MLE values for the operating model.
    tmp_total <- round(obj$env$data$total * sim_size,0)
    obj$env$data$total <- tmp_total
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
                         , silent = FALSE)
    
    sim_opt <- TMBhelper::fit_tmb( sim_obj,
                               loopnum = 1,
                               newtonsteps = 1,
                               quiet = TRUE,
                               getsd = FALSE) # where Obj is a compiled TMB object
    
    sim_obj$rep <- sim_obj$report()
    
    sim_out[i,] <- sim_opt$par
    sim_proj[i,,,] <- sim_obj$rep$proj_y
    
    print(paste("sys time", sys-Sys.time()))
  }
  # return(sim_out)
  return(list(sim_out = sim_out
              ,sim_i = sim_i
              ,sim_proj = sim_proj
              ,sim_obj = sim_obj
              ,sim_opt = sim_opt
              ))
}

