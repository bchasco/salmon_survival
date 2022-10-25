f_TMB_data_pars <- function(df = 0,
                            mesh = 0,
                            rangeL = 0,
                            rangeJ = 0,
                            AR_flags = 0,
                            by_pass = 0,
                            H_flag = 0,
                            b_flags = 0,
                            re_flags = 0){
  # library(RandomFields)
  library(raster)
  library(RANN)
  
  data <- list(n_i = length(df$s)
               ,a_i = df$a
               ,t_AR = AR_flags["t_AR"]
               ,l_AR = AR_flags["l_AR"]
               ,j_AR = AR_flags["j_AR"]
               ,mu_bypass = b_flags['mu_bypass']
               ,t_bypass = b_flags['t_bypass']
               ,j_bypass = b_flags['j_bypass']
               ,l_bypass = b_flags['l_bypass']
               ,jlt_bypass = b_flags['jlt_bypass']
               ,H_flag = H_flag
               ,j_flag = re_flags['j_flag']
               ,l_flag = re_flags['l_flag']
               ,t_flag = re_flags['t_flag']
               # ,z_jl_flag = z_jl_flag
               ,jlt_flag = re_flags['jlt_flag']
               ,n_t = max(df$y - min(df$y))+1
               ,x_s_jl = mesh$mesh_jl$idx$loc - 1
               ,s_i_jl = mesh$knots_xy_jl$cluster - 1
  #              ,management = cbind(s_i_jl_jm7_lp10 - 1,
  #                                  s_i_jl_jm7_lm10 - 1, 
  #                                  s_i_jl_jp7_lm10 - 1, 
  #                                  s_i_jl_jp7_lp10 - 1, 
  #                                  s_i_jl_jp7_l0 -1, 
  #                                  s_i_jl_jm7_l0 - 1, 
  #                                  s_i_jl_j0_lp10 -1 , 
  #                                  s_i_jl_j0_lm10 -1, 
  #                                  s_i_jl_j0_l0 -1)
               ,t_i = df$y - min(df$y)#Use all of the dataa
               ,l_i = df$l - rangeL['minL']
               ,j_i = df$j - rangeJ['minJ']
               ,total = df$nt
               ,surv = df$ns
               ,spde_jl = mesh$spde_jl
               ,n_s_jl = mesh$spde_jl$n_s
  )
  parameters <- list(
    mu = -4
    ,bypass = 0
    ,log_tau_jl2 = 0
    ,log_kappa_jl = 0
    ,ln_H_input_jl = c(0,0)
    # ,z_jl = rep(0,c(mesh$spde_jl$n_s)) #number of stations by the number of years
    ,z_jlt = array(0,c(mesh$spde_jl$n_s,
                       b_flags['jlt_bypass'] + 1,
                       data$n_t)) #num
    ,y_re = matrix(0, data$n_t, 
                   b_flags['t_bypass'] + 1)
    ,j_re = matrix(0, max(data$j_i)+1, 
                   b_flags['j_bypass'] + 1)
    ,l_re = matrix(0, max(data$l_i)+1, 
                   b_flags['l_bypass'] + 1)
    ,ln_sigl = 0
    ,ln_sigy = 0
    ,ln_sigj = 0
    ,f_phil = 0
    ,f_phiy = 0
    ,f_phij = 0
    ,f_psil = 0
    ,f_psiy = 0
    ,f_psij = 0
    ,f_psijl = 0
  )

  return(list(data = data,
              parameters = parameters))    
}  