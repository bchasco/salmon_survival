sqrt(1-(1/(1+exp(optList[[2]]$par["f_phij"])))^2) * exp(optList[[2]]$par["ln_sigj"])
sqrt(1-(1/(1+exp(optList[[2]]$par["f_phil"])))^2) * exp(optList[[2]]$par["ln_sigl"])
sqrt(1-(1/(1+exp(optList[[2]]$par["f_phiy"])))^2) * exp(optList[[2]]$par["ln_sigy"])
1/exp(optList[[2]]$par["log_tau2_jl"])
