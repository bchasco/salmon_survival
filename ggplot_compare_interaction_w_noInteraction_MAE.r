


par(mfrow=c(2,2))
for(i in 1:4){
  x1 <- simOut[!is.na(simOut[,1,1,1]),,i,1] #sim X year
  x2 <- simDat[!is.na(simOut[,1,1,1]),] #sim X yr
  
  plot(colMeans(abs(x1 - x2)), ylim=c(0,0.015), type="l")
  
  x1 <- simOut[!is.na(simOut[,1,1,1]),,i,2] #sim X year
  x2 <- simDat[!is.na(simOut[,1,1,1]),] #sim X yr
  
  lines(colMeans(abs(x1 - x2)))
  
}
