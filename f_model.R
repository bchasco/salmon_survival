f_model <- function(raw_file = raw_file,
                   rangeL = rangeL,
                   rangeJ = rangeJ,
                   AR_flags = AR_flags,
                   bypass_flags = bypass_flags,
                   H_flag = H_flag,
                   re_flags = re_flags,
                   version = version,
                   n_knots = n_knots,
                   random = random,
                   compare_AIC = compare_AIC,
                   getsd = getsd,
                   m_proj = proj,
                   save_to_file = save_to_file,
                   save_file = save_file,
                   DHARMa_sim = DHARMa_sim,
                   bias_sim = bias_sim,
                   bias_sim_n = bias_sim_n){
  
  #read the data and do any transformations
  df <- f_raw_data(file = raw_file,
                   rangeL = rangeL,
                   rangeJ = rangeJ)
  
  projections <- f_management_projection(df = df,
                                       m_proj = m_proj,
                                       rangeJ = rangeJ,
                                       rangeL = rangeL)
  #create the INLA mesh
  mesh <- 
    f_mesh(data = df,
                 n_knots = n_knots,
                 projections = projections)

  # #Create the data and parameters objects for TMB
  TMB_list <- 
  f_TMB_data_pars(df = df,
                              mesh = mesh,
                              by_pass = 1,
                              AR_flags = AR_flags,
                              H_flag = H_flag,
                              rangeL = rangeL,
                              rangeJ = rangeJ,
                              b_flags = bypass_flags,
                              re_flags = re_flags)

  #Create the necessary map
  map <- #NA
  f_map(TMB_list)

  #TMB object
  obj <-   #NA
  MakeADFun(TMB_list$data
                   ,silent = FALSE
                   , TMB_list$parameters
                   , random=random
                   , map = map
                   , DLL=paste0("spde_aniso_wt_",version)
                   , sdreport = FALSE)

  # print(obj$par)
  
  #keep track of the mesh
  # obj$mesh <- mesh

  #Estimate the parameters
  opt <- #NA
  TMBhelper::fit_tmb( obj,
                             loopnum = 1,
                             newtonsteps = 1,
                             quiet = TRUE,
                             getsd = getsd) # where Obj is a compiled TMB object

  #Save the report
  obj$rep <- obj$report()

  if(bias_sim){
    bias_sim <- f_bias_sim(obj = obj,
               TMB_list = TMB_list,
               map = map,
               bias_sim_n = bias_sim_n)
  }else{
    bias_sim <- NA
  }
  #This just keeps track of the model comparisons for 
  #choosing the best model
  # if(compare_AIC){
  #   load(file="AIC_comp.rData")
  #   tmp_list <- list(bypass_flags = bypass_flags,
  #                    AR_flags = AR_flags,
  #                    re_flags = re_flags,
  #                    AIC = round(opt$AIC,1))
  # 
  #   if(length(AIC_comp)==0){
  #     print(unlist(tmp_list))
  #     AIC_comp[[1]] <- unlist(tmp_list)
  #     save(file="AIC_comp.rData",AIC_comp)
  #   }
  #   if(length(AIC_comp)>=1){
  #     comps <- c('bypass_flags','re_flags','AR_flags')
  #     same <- FALSE
  #     for(ii in 1:length(AIC_comp)){
  #       if(sum(unlist(AIC_comp[[ii]])[names(AIC_comp[[ii]])!='AIC']!=unlist(tmp_list)[names(AIC_comp[[ii]])!='AIC'])==0){
  #         same <- TRUE
  #         print("You've already run this model.")
  #         tmp_df <- as.data.frame(AIC_comp, col.names=paste('model',1:length(AIC_comp)))
  #         print(tmp_df)
  #         print(paste("The delta AIC for this model is", round(abs(min(tmp_df['AIC',])-tmp_df['AIC',ii]),1)))
  #         print(paste("See model", ii))
  #         break;
  #       }
  #     }
  #     if(!same){
  #       print("not same")
  #       AIC_comp[[length(AIC_comp)+1]] <- unlist(tmp_list)
  #       save(file="AIC_comp.rData",AIC_comp)
  #       tmp_df <- as.data.frame(AIC_comp, 
  #                               col.names=paste('model',1:length(AIC_comp)))
  #       print(tmp_df)
  #     }
  #   }
  # }

  # d_res <- NA
  if(DHARMa_sim){
    
    sim <- replicate(100, {
      simdata <- obj$simulate()$surv
    })
    d_res <- createDHARMa(sim[,],
                          obj$env$data$surv,
                          fittedPredictedResponse = plogis(obj$rep$eta_i)*obj$env$data$total)
    KSi <- testUniformity(d_res)
  }
  
  if(save_to_file){
    fit <- list(
      obj = obj
      ,opt = opt
      ,mesh = mesh
      ,df = df
      ,TMB_list = TMB_list
      ,proj = projections)
    save(file=save_file,fit)
  }  
  return(list(
    obj = obj
    ,opt = opt
    ,mesh = mesh
    ,df = df
    ,TMB_list = TMB_list
    ,proj = projections
    # ,sim = sim
    # ,d_res = d_res
    # ,KSi = KSi
    ,bias_sim = bias_sim
))
}

