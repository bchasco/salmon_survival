par(mfrow=c(2,2))

load("modSD.rData")
load("modREP.rData")

library(ggplot2)

mar_y <- cbind(repList[[1]]$mar_y,repList[[2]]$mar_y,repList[[3]]$mar_y,repList[[4]]$mar_y)
mar_l <- cbind(repList[[1]]$mar_l,repList[[2]]$mar_l,repList[[3]]$mar_l,repList[[4]]$mar_l)
mar_j <- cbind(repList[[1]]$mar_j,repList[[2]]$mar_j,repList[[3]]$mar_j,repList[[4]]$mar_j)

matplot(mar_l)
matplot(mar_j)
matplot(mar_y)


mar_y <- cbind(sdList$sd[[1]]$mar_y,sdList$sd[[2]]$mar_y,sdList$sd[[3]]$mar_y,sdList$sd[[4]]$mar_y)
mar_l <- cbind(sdList$sd[[1]]$mar_l,sdList$sd[[2]]$mar_l,sdList$sd[[3]]$mar_l,sdList$sd[[4]]$mar_l)
mar_j <- cbind(sdList$sd[[1]]$mar_j,sdList$sd[[2]]$mar_j,sdList$sd[[3]]$mar_j,sdList$sd[[4]]$mar_j)

matplot(mar_l)
matplot(mar_j)
matplot(mar_y, type="l", lty=1)

par(mfrow=c(5,2))
for(yy in 11:19){
  mar_y <- cbind(sdList$est[[1]]$proj_surv[yy,],sdList$est[[2]]$proj_surv[yy,],sdList$est[[3]]$proj_surv[yy,],sdList$est[[4]]$proj_surv[yy,])
  matplot(mar_y, type="l", lty=1)
}
# mar_y <- cbind(sdList$sd[[1]]$proj_surv_ag,sdList$sd[[2]]$proj_surv_ag,sdList$sd[[3]]$proj_surv_ag,sdList$sd[[4]]$proj_surv_ag)
