f_TMB_data_pars <- function(df = 0, #pass the raw data
                            mesh = 0, #mesh
                            rangeL = 0, #range of the lengths
                            rangeJ = 0, #range of days
                            AR_flags = 0, #the rest are flags
                            by_pass = 0, 
                            H_flag = 0,
                            b_flags = 0,
                            re_flags = 0,
                            sim_size = 0,
                            proj_sim = 0,
                            proj_H = c(0,0)){
  
  library(raster)
  library(RANN)

  for(i in 1:dim(mesh$tmb_proj)[2]){
    mesh$tmb_proj[,i,'l'] <- mesh$tmb_proj[,i,'l'] - min(rangeL)
    mesh$tmb_proj[,i,'j'] <- mesh$tmb_proj[,i,'j'] - min(rangeJ)
  }
    
  data <- list(n_i = nrow(df)
               ,a_i = df$a
               ,t_AR = AR_flags["t_AR"]
               ,l_AR = AR_flags["l_AR"]
               ,j_AR = AR_flags["j_AR"]
               ,mu_bypass = b_flags['mu_bypass']
               ,t_bypass = b_flags['t_bypass']
               ,j_bypass = b_flags['j_bypass']
               ,l_bypass = b_flags['l_bypass']
               ,jl_bypass = b_flags['jl_bypass']
               ,jlt_bypass = b_flags['jlt_bypass']
               ,H_flag = H_flag
               ,j_flag = re_flags['j_flag']
               ,l_flag = re_flags['l_flag']
               ,t_flag = re_flags['t_flag']
               ,jl_flag = re_flags['jl_flag']
               ,jlt_flag = re_flags['jlt_flag']
               ,n_t = max(df$y - min(df$y))+1
               ,x_s_jl = mesh$mesh_jl$idx$loc - 1
               ,s_i_jl = mesh$knots_xy_jl$cluster - 1
               ,t_i = df$y - min(df$y)#Use all of the dataa
               ,l_i = df$l - rangeL['minL']
               ,j_i = df$j - rangeJ['minJ']
               ,total = df$nt
               ,surv = df$ns
               ,spde_jl = mesh$spde_jl
               ,n_s_jl = mesh$spde_jl$n_s
               ,m = mesh$tmb_proj
               ,m_imv = dim(mesh$tmb_proj)
               # ,sim_size = sim_size
               ,proj_sim = proj_sim
               ,proj_H = proj_H
  )
  parameters <- list(
    mu = -4
    ,bypass = 0
    ,log_tau_jl = 0
    ,log_tau_jl2 = 0
    ,log_kappa_jl = 0
    ,ln_H_input_jl = c(0,0)
    ,z_jl = array(0,c(mesh$spde_jl$n_s,
                      b_flags['jl_bypass'] + 1)) #number of stations by the number of years
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
    ,f_psijlt = 0
  )

  return(list(data = data,
              parameters = parameters))    
}  