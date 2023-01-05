f_plot_ROC <- function(fit = fit){
  library(pROC)
  library(ggplot2)
  
  #expand all of the predictions
  pred <- rep(plogis(fit$obj$rep$eta_i), fit$obj$env$data$total)
  
  #create two categories of surv and dead for each of the predictions
  obs <- NA
  for(i in 1:length(fit$obj$env$data$total)){
    obs <- c(obs,rep("surv",fit$obj$env$data$surv[i]),rep("dead",fit$obj$env$data$total[i]-fit$obj$env$data$surv[i]))
  }
  obs <- na.omit(obs)

  #create the roc   
  myRoc <- roc(obs,pred)

  print(myRoc)  
  png("f_plot_ROC.png", unit="in", height=5, width = 5, res = 600)
  rocobj <- plot.roc(obs, pred,
                     main = "Confidence intervals of specificity/sensitivity", 
                     percent = TRUE,
                     ci = TRUE, 
                     of = "se",                                # ci of sensitivity
                     specificities = seq(0, 100, 5),           # on a select set of specificities
                     ci.type="shape", 
                     ci.col="#1c61b6AA")      # plot the CI as a blue shape
  plot(rocobj)
  dev.off()
}