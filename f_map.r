f_map <- function(TMBlist = TMB_list){
  
  map <- list()

  if(TMBlist$data$mu_bypass == 0){
    map <- append(map,
                  list(bypass = as.factor(NA)))
  }
  
  if(TMBlist$data$t_flag==0){
    map <- append(map,
                  list(
                    ln_sigy = as.factor(NA)
                    ,f_phiy = as.factor(NA)
                    ,f_psiy = as.factor(NA)
                    ,y_re = as.factor(matrix(NA,nrow(TMBlist$parameters$y_re),
                                             ncol(TMBlist$parameters$y_re)
                  ))))
  }else{
    if(TMBlist$data$t_bypass==0){
      map <- append(map,
                    list(f_psiy = as.factor(NA)))
    }
    if(TMBlist$data$t_AR!=1){
      map <- append(map,
                    list(f_phiy = as.factor(NA)))
    }
  }

  
  if(TMBlist$data$j_flag==0){
    map <- append(map,
                  list(
                    ln_sigj = as.factor(NA)
                    ,f_phij = as.factor(NA)
                    ,f_psij = as.factor(NA)
                    ,j_re = as.factor(matrix(NA,nrow(TMBlist$parameters$j_re),
                                             ncol(TMBlist$parameters$j_re)
                    ))))
  }else{
    if(TMBlist$data$j_bypass==0){
      map <- append(map,
                    list(f_psij = as.factor(NA)))
    }
    if(TMBlist$data$j_AR!=1){
      map <- append(map,
                    list(f_phij = as.factor(NA)))
    }
  }

  if(TMBlist$data$l_flag==0){
    map <- append(map,
                  list(
                    ln_sigl = as.factor(NA)
                    ,f_phil = as.factor(NA)
                    ,f_psil = as.factor(NA)
                    ,l_re = as.factor(matrix(NA,
                                             nrow(TMBlist$parameters$l_re),
                                             ncol(TMBlist$parameters$l_re)
                    ))))
  }else{
    if(TMBlist$data$l_bypass==0){
      map <- append(map,
                    list(f_psil = as.factor(NA)))
    }
    if(TMBlist$data$l_AR!=1){
      map <- append(map,
                    list(f_phil = as.factor(NA)))
    }
  }

  #There's no AR1 model the jlt 
  if(TMBlist$data$jlt_flag == 0){
    jlt_dim <- dim(TMBlist$parameters$z_jlt)
    map <- append(map,
                  list(
                    log_tau_jl2 = as.factor(NA)
                    ,f_psijlt = as.factor(NA)
                    ,z_jlt = as.factor(array(NA,c(jlt_dim))) #num
                  ))
  }else{
    if(TMBlist$data$jlt_bypass == 0){
      map <- append(map,
                    list(
                      f_psijlt = as.factor(NA)
                    ))
    }
  }

  if(TMBlist$data$jl_flag == 0){
    jl_dim <- length(TMBlist$parameters$z_jl)
    map <- append(map,
                  list(
                    log_tau_jl = as.factor(NA)
                    ,f_psijl = as.factor(NA)
                    ,z_jl = as.factor(array(NA,c(jl_dim))) #num
                  ))
  }else{
    if(TMBlist$data$jl_bypass == 0){
      map <- append(map,
                    list(
                      f_psijl = as.factor(NA)
                    ))
    }
  }
  
  if(TMBlist$data$jlt_flag == 0 & TMBlist$data$jl_flag == 0){
    map <- append(map,
                  list(
                    log_kappa_jl = as.factor(NA)
                    ,ln_H_input_jl = as.factor(rep(NA,2))))
  }
  
  print(names(map))
  return(map)
  
}