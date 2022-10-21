simComp <- matrix(NA, nsim, 9)


for(i in 1:2){
  # simComp[i,] <- apply(plogis(simRep[[1]][[i]]$tmp[,]),2,median)/
  # apply(plogis(simRep[[2]][[i]]$tmp[,]),2,median)
  simComp[i,] <- plogis(simRep[[1]][[i]]$proj_surv_ag) /
    plogis(simRep[[2]][[i]]$proj_surv_ag)
}

par(mfrow=c(3,3))
for(i in 1:9){
  hist(simComp[1:nsim,i], 
       xlim=c(0.9,1.3), 
       las=1,
       ylab="Frequency",
       xlab="ratio",
       main=NULL)
  abline(v=1, col="red", lwd=2)
}
