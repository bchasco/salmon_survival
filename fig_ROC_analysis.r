library(pROC)
library(ggplot2)

load("sdList_v11.rData")
load("repList_v11.rData")
bm <- 3 #model with the lowest AIC


pred <- rep(plogis(repList[[bm]]$eta_i), data$total)

obs <- NA
for(i in 1:length(data$total)){
  obs <- c(obs,rep("surv",data$surv[i]),rep("dead",data$total[i]-data$surv[i]))
}
obs <- na.omit(obs)

myRoc <- roc(obs,pred)

png("fig_ROC_analysis.png", unit="in", height=5, width = 5, res = 600)
rocobj <- plot.roc(obs, pred,
                   main = "Confidence intervals of specificity/sensitivity", 
                   percent = TRUE,
                   ci = TRUE, 
                   of = "se",                                # ci of sensitivity
                   specificities = seq(0, 100, 5),           # on a select set of specificities
                   ci.type="shape", 
                   ci.col="#1c61b6AA")      # plot the CI as a blue shape
# plot(rocobj)
dev.off()
