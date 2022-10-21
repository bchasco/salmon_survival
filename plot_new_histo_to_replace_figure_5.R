xx <- matrix(NA,50,19)
# xx <- matrix(NA,50,19)
res <- matrix(NA,50,19)
res2 <- matrix(NA,50,19)
res3 <- matrix(NA,50,19)
opt$par[1]
rep <- org_obj$report()
for(i in c(1:50)){
  x1 <- sqrt(sum((simRep$GMRF[[i]]$j_re-rep$j_re)^2)/length(rep$j_re))
  x2 <- sqrt(sum((simRep$noGMRF[[i]]$j_re-rep$j_re)^2)/length(rep$j_re))
  # xx[i,] <- t(simRep$GMRF[[i]]$proj_surv_nox[9,]/
  #                   (simRep$GMRF[[i]]$proj_surv[9,]))
  for(ii in 1:19){
    # print(paste(c(1:9)[max(simRep$GMRF[[i]]$proj_surv[ii,])==simRep$GMRF[[i]]$proj_surv[ii,]]==
    #               c(1:9)[max(simRep$noGMRF[[i]]$proj_surv[ii,])==simRep$noGMRF[[i]]$proj_surv[ii,]]))
    xx[i,ii] <- c(1:9)[max(simRep$GMRF[[i]]$proj_surv[ii,])==simRep$GMRF[[i]]$proj_surv[ii,]]
    res[i,ii] <- plogis(simRep$GMRF[[i]]$proj_surv[ii,xx[i,ii]])/plogis(simRep$GMRF[[i]]$proj_surv[ii,1])
    res2[i,ii] <- plogis(simRep$GMRF[[i]]$proj_surv[ii,1])/plogis(simRep$GMRF[[i]]$proj_surv[ii,9])
    res3[i,ii] <- plogis(simRep$noGMRF[[i]]$proj_surv[ii,1])/plogis(simRep$noGMRF[[i]]$proj_surv[ii,9])
  }
}

par(mfrow=c(3,1))
hx <- hist(na.omit(res2[]), breaks=50, col="red", ylim=c(0,100), xlab=""); 
hist(na.omit(res3[xx==1]),add=TRUE, breaks=hx$breaks, col="blue", xlab="Percent change in survival"); 

hist(na.omit(xx[xx!=1]), xlab="Management scenario"); 
hist(na.omit(res[xx!=1]), xlab="Percent change in survival"); 

