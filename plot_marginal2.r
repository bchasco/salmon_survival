load("repList_v11.rData")
load("sdList_v11.rData")
par(mfrow=c(2,2))
   plot((plogis(sdList$est[[3]]$y_est)-plogis(sdList$est[[1]]$y_est))/(plogis(sdList$est[[1]]$mar_y))*100, las = 1, type="l")
   plot(plogis(sdList$est[[1]]$mar_y), type="l", ylim=c(0,0.08))
 lines(plogis(sdList$est[[3]]$mar_y), col="red", lty=2)
 lines(plogis(sdList$est[[1]]$mar_y+1.28*sdList$sd[[1]]$mar_y))
 lines(plogis(sdList$est[[1]]$mar_y-1.28*sdList$sd[[1]]$mar_y))
 lines(plogis(sdList$est[[3]]$mar_y+1.28*sdList$sd[[3]]$mar_y), col="red")
 lines(plogis(sdList$est[[3]]$mar_y-1.28*sdList$sd[[3]]$mar_y), col="red")
 
 
          plogis(sdList$est[[3]]$mar_y[i]), 
          plogis(sdList$est[[1]]$mar_y[i]-1.28*sdList$sd[[1]]$mar_y[i]),
          plogis(sdList$est[[3]]$mar_y[i]),
          col="lightgrey")
 
 lm(plogis(sdList$est[[1]]$mar_y)~-1+plogis(sdList$est[[3]]$mar_y))
 
 for(i in 1:length(sdList$est[[3]]$mar_y)){
   segments(plogis(sdList$est[[1]]$mar_y[i]),
            plogis(sdList$est[[3]]$mar_y[i]-1.28*sdList$sd[[3]]$mar_y[i]), 
            plogis(sdList$est[[1]]$mar_y[i]),
            plogis(sdList$est[[3]]$mar_y[i]+1.28*sdList$sd[[3]]$mar_y[i]),
            col="lightgrey")
 }
 segments(0,0,1,1)

 
 plot(plogis(sdList$est[[1]]$mar_j))
 lines(plogis(sdList$est[[3]]$mar_j), pch=16, cex=1.4, xlim=c(0,0.02),ylim=c(0,0.02))
 for(i in 1:length(sdList$est[[3]]$mar_j)){
   segments(plogis(sdList$est[[1]]$mar_j[i]+1.28*sdList$sd[[1]]$mar_j[i]),
            plogis(sdList$est[[3]]$mar_j[i]), 
            plogis(sdList$est[[1]]$mar_j[i]-1.28*sdList$sd[[1]]$mar_j[i]),
            plogis(sdList$est[[3]]$mar_j[i]),
            col="lightgrey")
 }
 
 for(i in 1:length(sdList$est[[3]]$mar_j)){
   segments(plogis(sdList$est[[1]]$mar_j[i]),
            plogis(sdList$est[[3]]$mar_j[i]-1.28*sdList$sd[[3]]$mar_j[i]), 
            plogis(sdList$est[[1]]$mar_j[i]),
            plogis(sdList$est[[3]]$mar_j[i]+1.28*sdList$sd[[3]]$mar_j[i]),
            col="lightgrey")
 }
 segments(0,0,1,1)
 points(plogis(sdList$est[[1]]$mar_j),plogis(sdList$est[[3]]$mar_j), pch=16, cex=1.4, xlim=c(0,0.02),ylim=c(0,0.02))
 lm(plogis(sdList$est[[1]]$mar_j)~-1+plogis(sdList$est[[3]]$mar_j))
 
 
 plot(plogis(sdList$est[[1]]$mar_l))
 lines(plogis(sdList$est[[3]]$mar_l), pch=16, cex=1.4, xlim=c(0,0.02),ylim=c(0,0.02))
 for(i in 1:length(sdList$est[[3]]$mar_l)){
   segments(plogis(sdList$est[[1]]$mar_l[i]+1.28*sdList$sd[[1]]$mar_l[i]),
            plogis(sdList$est[[3]]$mar_l[i]), 
            plogis(sdList$est[[1]]$mar_l[i]-1.28*sdList$sd[[1]]$mar_l[i]),
            plogis(sdList$est[[3]]$mar_l[i]),
            col="lightgrey")
 }
 
 for(i in 1:length(sdList$est[[3]]$mar_l)){
   segments(plogis(sdList$est[[1]]$mar_l[i]),
            plogis(sdList$est[[3]]$mar_l[i]-1.28*sdList$sd[[3]]$mar_l[i]), 
            plogis(sdList$est[[1]]$mar_l[i]),
            plogis(sdList$est[[3]]$mar_l[i]+1.28*sdList$sd[[3]]$mar_l[i]),
            col="lightgrey")
 }
 segments(0,0,1,1)
 points(plogis(sdList$est[[1]]$mar_l),plogis(sdList$est[[3]]$mar_l), pch=16, cex=1.4, xlim=c(0,0.02),ylim=c(0,0.02))
 lm(plogis(sdList$est[[1]]$mar_l)~-1+plogis(sdList$est[[3]]$mar_l))
 