par(mfrow=c(3,4))
for(i in 1:length(fit$opt$par)){
  hist(mySim[,i], main = names(fit$opt$par[i]))
  abline(v = fit$opt$par[i], col="red", lwd=2)
}