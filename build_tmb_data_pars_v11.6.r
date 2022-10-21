# library(RandomFields)
library(raster)
library(RANN)

data <- list(n_i = n_i
             ,b_i = df2$a
             ,n_t = length(unique(df2$y))
             ,t_i = df2$y - min(df2$y)#Use all of the dataa
             ,l_i = df2$len - minL
             ,j_i = df2$j - minJ
             ,total = df2$ns
             ,surv = df2$s
)


data$spde_jl <- spde_jl
data$n_s_jl <- data$spde_jl$n_s


parameters <- list(
  mu = -4
  ,bypass = 0
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
)
