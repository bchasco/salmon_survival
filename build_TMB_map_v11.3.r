map <- list()

if(!data$z_jl_flag){
  map <- append(map,
                list(
                  log_tau_jl = as.factor(NA)
                  ,z_jl = as.factor(rep(NA,c(data$spde_jl$n_s)))
                )
  )
}


if(!data$z_jl_flag & !data$z_jlt_flag){
  map <- append(map,
                list(
                  ln_H_input_jl = as.factor(rep(NA,2))
                  ,log_kappa_jl = as.factor(NA)
                )
  )
}


if(!data$z_jlt_flag){
  map <- append(map,
                list(
                  log_tau_jl2 = as.factor(NA)
                  # ,f_phijlt = as.factor(NA)
                  ,z_jlt = as.factor(matrix(NA,data$spde_jl$n_s, nrow(parameters$y_re)))
                )
  )
}

# if(!data$z_jt_flag){
#   map <- append(map,
#                 list(
#                   f_phijt_2 = as.factor(NA)
#                   ,ln_sigjt = as.factor(NA)
#                   ,z_jt = as.factor(matrix(NA,dim(parameters$z_jt)[1], dim(parameters$z_jt)[2]))
#                 )
#   )
# }

if(!data$t_flag){
  map <- append(map,
                list(
                  ln_sigy = as.factor(NA)
                  ,f_phiy = as.factor(NA)
                  ,y_re = as.factor(rep(NA,data$n_t))
                )
  )
}

if(!data$j_flag){
  map <- append(map,
                list(
                  ln_sigj = as.factor(NA)
                  ,f_phij = as.factor(NA)
                  ,j_re = as.factor(rep(NA,length(parameters$j_re)))
                )
  )
}

if(!data$l_flag){
  map <- append(map,
                list(
                  ln_sigl = as.factor(NA)
                  ,f_phil = as.factor(NA)
                  ,l_re = as.factor(rep(NA,length(parameters$l_re)))
                )
  )
}

