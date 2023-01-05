f_bias_sim <- function(df = df,
                       opt = opt, 
                       obj = obj,
                       TMB_list = TMB_list,
                       map = map,
                       lwr = lwr,
                       upr = upr,
                       bias_sim_n = bias_sim_n,
                       sim_size = sim_size){
  
  
  sim_out <- matrix(NA,bias_sim_n,length(obj$par))
  sim_proj <- array(NA,c(bias_sim_n,dim(obj$rep$proj_y)))
  sim_y <- array(NA,c(bias_sim_n,length(unique(obj$env$data$t_i)),3))
  
  obj$env$data$sim_size <- sim_size
  obj$env$data$total <- round(obj$env$data$total * sim_size,0)
  # fit$obj$env$parList()$mu <- -4.5
    for(i in 1:bias_sim_n){
    sys <- Sys.time()
    #Set the year, day, and length random effects equal to their MLE values for the operating model.
    # obj$env$last.par['mu'] <- (-4.5)
    sim_i <- obj$simulate(complete=TRUE) #simulate data from the original model
    print(paste("sim ", i , "out of ", bias_sim_n))
    random <- obj$env$.random
    sim_obj <- MakeADFun(data=sim_i #create a new TMB model with the simulated data
                         , parameters=TMB_list$parameters
                         , random=random
                         , map = map
                         , DLL=paste0("spde_aniso_wt_","v11_6")
                         # , sdreport = FALSE
                         , checkParameterOrder = TRUE
                         , silent = TRUE)
    print(paste("Simulated total,", sum(sim_obj$env$data$total)))  
    print(paste("Observed survivors,",sum(obj$env$data$surv[obj$env$data$a_i==0]),sum(obj$env$data$surv[obj$env$data$a_i==1])))  
    print(paste("Simulated survivors,",sum(sim_obj$env$data$surv[sim_obj$env$data$a_i==0]),sum(sim_obj$env$data$surv[sim_obj$env$data$a_i==1])))  
    
    sim_opt <- TMBhelper::fit_tmb( sim_obj,
                               loopnum = 1,
                               newtonsteps = 1,
                               startpar = opt$par,
                               lower = rep(-6.5,length(sim_obj$par)),
                               upper = rep(3,length(sim_obj$par)),
                               quiet = TRUE,
                               getsd = FALSE) # where Obj is a compiled TMB object
    
    sim_obj$rep <- sim_obj$report()
    
    sim_out[i,] <- sim_opt$par
    sim_proj[i,,,] <- sim_obj$rep$proj_y
    
    sim_df <- data.frame(y = sim_obj$env$data$t_i,
                         a = sim_obj$env$data$a_i,
                         nt = sim_obj$env$data$total,
                         ns = sim_obj$env$data$surv)
    
    tmp <- sim_df %>% 
      group_by(a = a, y = y+1998) %>% 
      summarise(p = round(sum(ns)/sum(nt),5)) %>% 
      pivot_wider(names_from = a, values_from = c(p), values_fill = 0) %>%
      arrange(y)
    
    # print(t(tmp))
    # print(i)
    # print(sim_y[i,,])
    # 
    # sim_y[i,,] <- unclass(tmp)
    
    print("simulation estimate")
    print(sim_opt$par)
    print("'true' estimate")
    print(opt$par)
    print(paste("sys time", round(Sys.time()-sys,2)))
  }
  # return(sim_out)
  return(list(sim_out = sim_out
              ,sim_i = sim_i
              ,sim_y = sim_y
              ,sim_proj = sim_proj
              ,sim_obj = sim_obj
              ,sim_opt = sim_opt
              ))
}

