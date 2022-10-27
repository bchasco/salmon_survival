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
                   getsd = getsd){
  
  #read the data and do any transformations
  df <- f_raw_data(file = raw_file,
                   rangeL = rangeL,
                   rangeJ = rangeJ)
  
  #create the INLA mesh
  mesh <- f_mesh(data = df,
                 n_knots = n_knots)
  
  #Create the data and parameters objects for TMB
  TMB_list <- f_TMB_data_pars(df = df,
                              mesh = mesh,
                              by_pass = 1,
                              AR_flags = AR_flags,
                              H_flag = H_flag,
                              rangeL = rangeL,
                              rangeJ = rangeJ,
                              b_flags = bypass_flags,
                              re_flags = re_flags)
  
  #Create the necessary map
  map <- f_map(TMB_list)
  
  #TMB object
  obj <- MakeADFun(TMB_list$data
                   ,silent = FALSE
                   , TMB_list$parameters
                   , random=random
                   , map = map
                   , DLL=paste0("spde_aniso_wt_",version)
                   , sdreport = FALSE)

  print(obj$par)
  
  #keep track of the mesh
  obj$mesh <- mesh

  #Estimate the parameters
  opt <- TMBhelper::fit_tmb( obj,
                             loopnum = 1,
                             newtonsteps = 1,
                             quiet = TRUE,
                             getsd = getsd) # where Obj is a compiled TMB object

  #Save the report
  obj$rep <- obj$report()

  #This just keeps track of the model comparisons for 
  #choosing the best model
  if(compare_AIC){
    load(file="AIC_comp.rData")
    tmp_list <- list(bypass_flags = bypass_flags,
                     AR_flags = AR_flags,
                     re_flags = re_flags,
                     AIC = round(opt$AIC,1))

    if(length(AIC_comp)==0){
      print(unlist(tmp_list))
      AIC_comp[[1]] <- unlist(tmp_list)
      save(file="AIC_comp.rData",AIC_comp)
    }
    if(length(AIC_comp)>=1){
      comps <- c('bypass_flags','re_flags','AR_flags')
      same <- FALSE
      for(ii in 1:length(AIC_comp)){
        if(sum(unlist(AIC_comp[[ii]])[names(AIC_comp[[ii]])!='AIC']!=unlist(tmp_list)[names(AIC_comp[[ii]])!='AIC'])==0){
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
  
  return(list(
    obj = obj
    ,opt = opt
    ,mesh = mesh
    ,df = df
    ,TMB_list = TMB_list))
}

