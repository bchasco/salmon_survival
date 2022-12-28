f_model <- function(raw_file = raw_file,
                   rangeL = NA,
                   rangeJ = NA,
                   AR_flags = NA,
                   bypass_flags = NA,
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
                   bias_sim_n = bias_sim_n,
                   sim_size = sim_size,
                   proj_sim = proj_sim,
                   proj_H = proj_H){
  
  #read the data and do any transformations
  df <- f_raw_data(file = raw_file,
                   rangeL = rangeL,
                   rangeJ = rangeJ)
  
  print("test")  
  projections <- f_management_projection(df = df,
                                       m_proj = m_proj,
                                       rangeJ = rangeJ,
                                       rangeL = rangeL)
  
  #create the INLA mesh
  mesh <- f_mesh(data = df,
                 n_knots = n_knots,
                 projections = projections)

  
  # #Create the data and parameters objects for TMB
  TMB_list <- f_TMB_data_pars(df = df,
                              mesh = mesh,
                              by_pass = 1,
                              AR_flags = AR_flags,
                              H_flag = H_flag,
                              rangeL = rangeL,
                              rangeJ = rangeJ,
                              b_flags = bypass_flags,
                              re_flags = re_flags,
                              sim_size = sim_size)
  #Create the necessary map
  # print(dim(df))
  # return(df)
  map <- #NA
  f_map(TMB_list)

  map_names <- names(map)
  map_names <- map_names[!(map_names%in%c(grep("_re",map_names),grep("z_",map_names)))]

  #fixed effect parameter length  
  par_names <- names(TMB_list$parameters)
  par_names <- par_names[!(par_names%in%c(grep("_re",par_names),grep("z_",par_names)))]
  par_names <- par_names[!(par_names%in%map_names)]
  
  print("par_names")
  print(par_names[order(par_names)])
  print("map_names")
  print(map_names[order(map_names)])
  
  upr <- rep(4,length(par_names))
  lwr <- rep(-6,length(par_names))

  #TMB object
  obj <-   #NA
  MakeADFun(TMB_list$data
                   ,silent = FALSE
                   , TMB_list$parameters
                   , random=random
                   , map = map
                    ,lower = lwr
                    ,upper = upr
                   , DLL=paste0("spde_aniso_wt_",version)
                   , sdreport = FALSE)
  print(obj$par)
  # print(names(map))
  
  #keep track of the mesh
  obj$mesh <- mesh
  #Estimate the parameters
  if(bias_sim) getsd <- FALSE
  opt <- #NA
  TMBhelper::fit_tmb( obj,
                     loopnum = 1,
                     newtonsteps = 1,
                    lower = rep(-6.5,length(obj$par)),
                    upper = rep(3,length(obj$par)),
                     quiet = TRUE,
                     getsd = getsd) # where Obj is a compiled TMB object

  #Save the report
  obj$rep <- obj$report()

  if(bias_sim){
    bias_sim <- f_bias_sim(opt = opt,
                           obj = obj,
               TMB_list = TMB_list,
               map = map,
               lwr = lwr,
               upr = upr,
               bias_sim_n = bias_sim_n,
               sim_size = sim_size)
  }else{
    bias_sim <- NA
  }
  #This just keeps track of the model comparisons for
  #choosing the best model
  if(compare_AIC){
    load(file="AIC_comp.rData")
    tmp_list <- list(bypass_flags = bypass_flags,
                     AR_flags = AR_flags,
                     re_flags = re_flags,
                     n_knots = n_knots,
                     AIC = round(opt$AIC,1))

    if(length(AIC_comp)==0){
      print(unlist(tmp_list))
      AIC_comp[[1]] <- unlist(tmp_list)
      save(file="AIC_comp.rData",AIC_comp)
    }
    if(length(AIC_comp)>=1){
      comps <- c('bypass_flags','re_flags','AR_flags','n_knots')
      same <- FALSE
      for(ii in 1:length(AIC_comp)){
        old_mod <- unlist(AIC_comp[[ii]])[names(AIC_comp[[ii]])!='AIC']
        new_mod <- unlist(tmp_list)[names(AIC_comp[[ii]])!='AIC'] 
        if(sum(old_mod!=new_mod)==0){
          same <- TRUE
          print("You've already run this model.")
          tmp_df <- as.data.frame(AIC_comp, col.names=paste('model',1:length(AIC_comp)))
          print(tmp_df)
          print(paste("The delta AIC for this model is", round(abs(min(tmp_df['AIC',])-tmp_df['AIC',ii]),1)))
          print(paste("See model", ii))
          break;
        }
      }
      if(!same){
        print("not same")
        AIC_comp[[length(AIC_comp)+1]] <- unlist(tmp_list)
        save(file="AIC_comp.rData",AIC_comp)
        tmp_df <- as.data.frame(AIC_comp,
                                col.names=paste('model',1:length(AIC_comp)))
        print(tmp_df)
      }
    }
  }

  d_res <- NA
  if(DHARMa_sim){

    dH_sim <- replicate(250, {
      simdata <- obj$simulate()$surv
    })
    d_res <- createDHARMa(dH_sim[,],
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
      ,bias_sim = bias_sim
      ,proj = projections)
    save(file=save_file,fit)
  }
  return(list(
    obj = obj
    ,opt = opt
    ,mesh = mesh
    ,df = df
    ,map = map
    ,TMB_list = TMB_list
    ,proj = projections
    ,bias_sim = bias_sim
))
}

