myMap <- list()
if(length(b)==1){
  myMap_b <- list(
    frho_y2  = as.factor(NA)
    ,frho_q  = as.factor(NA)
  )
  myMap <- append(myMap,
                  myMap_b)
}

if(re_y==0)
{
  myMap_y <- list(
    eps_y = as.factor(matrix(NA,1,data$n_t))
    ,frho_y  = as.factor(NA)
    ,ln_sigma_y  = as.factor(NA)
  )
  myMap <- append(myMap,
                  myMap_y)
}
if(s==0){
  myMap_s <- list(
    log_tau_O = as.factor(NA)
    ,log_kappa = as.factor(NA)
    ,Omega_input  = as.factor(rep(NA,n_x))
  )
  myMap <- append(myMap,
                  myMap_s)
}
if(st==0){
  myMap_st <- list(
    log_tau_E = as.factor(NA)
    ,frho = as.factor(NA)
    ,Epsilon_input  = as.factor(matrix(NA,n_x, data$n_t))
  )
  myMap <- append(myMap,
                  myMap_st)
  
}

mod <- paste0("st_",st,
              "_s_",s,
              "_rey_",re_y,
              "_dBin_",dBin,
              "_lBin_",lBin,
              "_b_",paste0(b,collapse = ""),
              "_maxBindl_",paste0(max_d,max_l),
              "_Datadlim_",paste0(min_data_d,max_data_d,collapse=""),
              "_Datallim_",paste0(min_data_l,max_data_l,collapse="")
)
