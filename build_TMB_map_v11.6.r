map <- list()

if(!data$t_flag){
  if(data$t_bypass){
    map <- append(map,
                  list(
                    ln_sigy = as.factor(NA)
                    ,f_phiy = as.factor(NA)
                    ,f_psiy = as.factor(NA)
                    ,y_re = as.factor(matrix(NA,2,data$n_t))
                  )
    )
  }else{
    map <- append(map,
                  list(
                    ln_sigy = as.factor(NA)
                    ,f_psiy = as.factor(NA)
                    ,f_phiy = as.factor(NA)
                    ,y_re = as.factor(rep(NA,data$n_t))
                  )
    )
    parameters$y_re <- matrix(0,data$n_t,1)
  }
}
if(!data$t_bypass){
  map <- append(map,
                list(f_psiy = as.factor(NA))
  )
  parameters$y_re <- matrix(0,data$n_t,1)
}



if(!data$j_flag){
  map <- append(map,
                list(
                  ln_sigj = as.factor(NA)
                  ,f_phij = as.factor(NA)
                  ,f_psij = as.factor(NA)
                  ))
  if(data$j_bypass){
    map <- append(map,
                  list(
                    j_re = as.factor(matrix(NA,2,max(data$j_i)+1))
                  ))
  }else{
    map <- append(map,
                  list(
                    j_re = as.factor(rep(NA,max(data$j_i)+1))
                  ))
    parameters$j_re <- matrix(0,max(data$j_i)+1,1)
  }
}
if(!data$j_bypass){
  map <- append(map,
                list(f_psij = as.factor(NA))
  )
  parameters$j_re <- matrix(0,max(data$j_i)+1,1)
}

if(!data$l_flag){
  map <- append(map,
                list(
                  ln_sigl = as.factor(NA)
                  ,f_phil = as.factor(NA)
                  ,f_psil = as.factor(NA)
                ))
  if(data$l_bypass){
    map <- append(map,
                  list(
                    l_re = as.factor(matrix(NA,2,max(data$l_i)+1))
                  ))
  }else{
    map <- append(map,
                  list(
                    l_re = as.factor(rep(NA,max(data$l_i)+1))
                  ))
    parameters$l_re <- matrix(0,max(data$l_i)+1,1)
  }
}
if(!data$l_bypass){
  map <- append(map,
                list(f_psil = as.factor(NA))
  )
  parameters$l_re <- matrix(0,max(data$l_i)+1,1)
}
