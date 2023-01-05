f_plot_DHARMa <- function(fit = fit,
                          nsample = 250){

  sim <- replicate(nsample, {
    simdata <- fit$obj$simulate(complete = TRUE)
    simdata <- simdata$surv
  })
  
  d_res <- createDHARMa(sim[,],
                        fit$obj$env$data$surv,
                        integerResponse = TRUE,
                        fittedPredictedResponse = plogis(fit$obj$rep$eta_i)*fit$obj$env$data$total)
  KSi <- testUniformity(d_res, plot = TRUE)
  # plotResiduals(d_res)
  # testDispersion(d_res)
  
}
