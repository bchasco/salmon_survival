# library(RandomFields)
library(raster)
library(RANN)

data <- list(n_i = n_i
             ,a_i = df2$a
             ,n_t = length(unique(df2$y))
             ,x_s_jl = mesh_jl$idx$loc - 1
             ,s_i_jl = s_i_jl - 1
             ,management = cbind(s_i_jl_jm7_lp10 - 1,
                                 s_i_jl_jm7_lm10 - 1, 
                                 s_i_jl_jp7_lm10 - 1, 
                                 s_i_jl_jp7_lp10 - 1, 
                                 s_i_jl_jp7_l0 -1, 
                                 s_i_jl_jm7_l0 - 1, 
                                 s_i_jl_j0_lp10 -1 , 
                                 s_i_jl_j0_lm10 -1, 
                                 s_i_jl_j0_l0 -1)
             ,t_i = df2$y - min(df2$y)#Use all of the dataa
             ,l_i = df2$len - minL
             ,j_i = df2$j - minJ
             ,env = matrix(0,1,1)
             ,total = df2$ns
             ,surv = df2$s
)


data$spde_jl <- spde_jl
data$n_s_jl <- data$spde_jl$n_s

data$j_im7 <- data$j_i - 7
data$j_im7[data$j_im7<0] <- 0
data$j_ip7 <- data$j_i + 7
data$j_ip7[data$j_ip7>max(data$j_i)] <- max(data$j_i)

data$j_shift <- cbind(data$j_im7,data$j_im7,data$j_ip7,data$j_ip7, data$j_ip7,data$j_im7, data$j_i, data$j_i, data$j_i)

data$l_im10 <- data$l_i - 4
data$l_im10[data$l_im10<0] <- 0
data$l_ip10 <- data$l_i + 4
data$l_ip10[data$l_ip10>max(data$l_i)] <- max(data$l_i)

data$l_shift <- cbind(data$l_ip10,data$l_im10,data$l_im10,data$l_ip10, data$l_i,data$l_i, data$l_ip10, data$l_im10, data$l_i)

parameters <- list(
  mu = -4
  ,bypass = 0
  ,log_tau_jt = 0
  ,log_tau_jl = 0
  ,log_tau_jl2 = 0
  ,log_kappa_jl = 0
  ,ln_H_input_jl = c(0,0)
  ,z_jl = rep(0,c(data$spde_jl$n_s)) #number of stations by the number of years
  ,z_jlt = array(0,c(data$spde_jl$n_s,2,data$n_t)) #num
  ,z_jt = matrix(0,max(data$j_i)+1,max(data$t_i)+1) #num
  ,y_re = matrix(0, data$n_t,2)
  ,j_re = matrix(0, max(data$j_i)+1,2)
  ,l_re = matrix(0, max(data$l_i)+1,2)
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

data <- append(data,
               list(
                 sim_lre = rep(0,length(parameters$l_re))
                 ,sim_lre_sd = rep(0,length(parameters$l_re))
                 ,sim_jre = rep(0,length(parameters$j_re))
                 ,sim_jre_sd = rep(0,length(parameters$j_re))))
  
  