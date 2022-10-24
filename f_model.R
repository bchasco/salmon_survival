f_model <- function(raw_file = raw_file,
                   rangeL = rangeL,
                   rangeJ = rangeJ,
                   bypass_flags = bypass_flags,
                   H_flag = H_flag,
                   re_flags = re_flags,
                   version = version,
                   n_knots = n_knots,
                   random = random){
  
  df <- f_raw_data(file = raw_file,
                   rangeL = rangeL,
                   rangeJ = rangeJ)
  
  mesh <- f_mesh(data = df,
                 n_knots = n_knots)
  
  TMB_list <- f_TMB_data_pars(df = df,
                              mesh = mesh,
                              by_pass = 1,
                              H_flag = H_flag,
                              rangeL = rangeL,
                              rangeJ = rangeJ,
                              b_flags = bypass_flags,
                              re_flags = re_flags)
  
  print(names(TMB_list))
  map <- f_map(TMB_list)
  
  obj <- MakeADFun(TMB_list$data
                   ,silent = FALSE
                   , TMB_list$parameters
                   , random=random
                   , map = map
                   , DLL=paste0("spde_aniso_wt_",version)
                   , sdreport = FALSE)
  
  obj$mesh <- mesh
  
  opt <- TMBhelper::fit_tmb( obj,
                             loopnum = 1,
                             newtonsteps = 1,
                             quiet = TRUE,
                             getsd = FALSE) # where Obj is a compiled TMB object
  
  obj$rep <- obj$report()
  
  
}

