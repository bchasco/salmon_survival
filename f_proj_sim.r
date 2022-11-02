sig <- 1/exp(fit$obj$env$last.par['log_tau_jl2'])
sig <- sig*c(0.5,1,2)
proj_tau <- log(1/sig)
f_tau_simulation <- function(fit,
                             proj_tau = proj_tau,
                             proj_sim_n = 20,
                             seed = 100){
  
  proj_sim <- array(NA, c(proj_sim_n,3,c(dim(fit$obj$rep$proj_y))))
  jcnt <- 1
  for(j in proj_tau){
    for(i in 1:proj_sim_n){
      fit$obj$env$data$proj_sim <- 1
      fit$obj$env$last.par[which(names(fit$obj$env$last.par) == 'log_tau_jl2')] <- j
      proj_sim[i,jcnt,,,] <- fit$obj$simulate(complete=TRUE)$proj_y
    }
    jcnt <- jcnt + 1
  }
  # 
  # tmp <- na.omit(c(proj_sim[,,1:2,3,])) 
  # hist(tmp, breaks=50)
  # print(sum(tmp<0)/length(tmp))
  # 
  # abline(v=0,lwd=2,col="red")
  return(list(proj_sim = proj_sim))
}
tau_proj <- f_tau_simulation(fit = fit)

