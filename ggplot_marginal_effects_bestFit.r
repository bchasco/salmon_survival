f_ggplot_marginal_effects <- function(fit){

  library(ggplot2)
  library(ggpubr)
  
  sdEst <- sdList$est[[bm]]
  sdSd <- sdList$sd[[bm]]
  mar_re <- c("mar_l","mar_j","mar_y")
  
  l_obs <- aggregate(list(wt_s = data$surv, t=data$total),by=list(l=data$l_i),sum)      
  j_obs <- aggregate(list(wt_s = data$surv, t=data$total),by=list(j=data$j_i),sum)      
  y_obs <- aggregate(list(wt_s = data$surv, t=data$total),by=list(y=data$t_i),sum)      
  
  myDf <- data.frame(var=c(rep(mar_re[1],length(sdEst$mar_l)),rep(mar_re[2],length(sdEst$mar_j)),rep(mar_re[3],length(sdEst$mar_y))),
                     x = c(1:length(parameters$l_re)+minL-1,1:length(parameters$j_re)+min(x$julian)-1,min(x$year):max(x$year)),
                     cat = rep("mar",length(c(sdEst$mar_l,sdEst$mar_j,sdEst$mar_y))),
                     val = c(sdEst$mar_l,sdEst$mar_j,sdEst$mar_y),
                     sd=c(sdSd$mar_l,sdSd$mar_j,sdSd$mar_y),
                     obs = c(l_obs$wt_s/l_obs$t, j_obs$wt_s/j_obs$t, y_obs$wt_s/y_obs$t ))
  myDf <- rbind(myDf,
                data.frame(var=c(rep(mar_re[1],length(sdEst$mar_l)),rep(mar_re[2],length(sdEst$mar_j)),rep(mar_re[3],length(sdEst$mar_y))),
                           x = c(1:length(parameters$l_re)+minL-1,1:length(parameters$j_re)+min(x$julian)-1,min(x$year):max(x$year)),
                           cat = rep("est",length(c(sdEst$l_est,sdEst$j_est,sdEst$proj_surv[,9]))),
                           val = c(sdEst$l_est,sdEst$j_est,sdEst$y_est),
                           sd=c(sdSd$l_est,sdSd$j_est,sdSd$y_est),
                           obs = c(l_obs$wt_s/l_obs$t, j_obs$wt_s/j_obs$t, y_obs$wt_s/y_obs$t )))
  
  myDf <- rbind(myDf,
                data.frame(var=c(rep(mar_re[1],length(sdEst$mar_l)),rep(mar_re[2],length(sdEst$mar_j)),rep(mar_re[3],length(sdEst$mar_y))),
                           x = c(1:length(parameters$l_re)+minL-1,1:length(parameters$j_re)+min(x$julian)-1,min(x$year):max(x$year)),
                           cat = rep("re",length(c(sdEst$l_re,sdEst$j_re,sdEst$y_re))),
                           val = c(sdEst$l_re,sdEst$j_re,sdEst$y_re),
                           sd=c(sdSd$l_re,sdSd$j_re,sdSd$y_re),
                           obs = c(l_obs$wt_s/l_obs$t, j_obs$wt_s/j_obs$t, y_obs$wt_s/y_obs$t )))
  
  
  g <- list(re=list(), est=list())
  myCat <- c("re","est")
  xlab <- c("Length (mm)", "Calendar day", "Migration year")
  for(i in 1:3){
    for(j in 1:2){
      
      g[[j]][[i]] <- ggplot(myDf[myDf$var==mar_re[i] & myDf$cat==myCat[j],], aes(x=x, y=val)) +
        theme_bw() + 
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
        ylim(0,0.038)
      if(j==1){
        g[[j]][[i]] <- g[[j]][[i]] +
          geom_line(aes(x=x, y=plogis(val + repList[[bm]]$mu))) +
          geom_ribbon(aes(ymin=plogis(val  + repList[[bm]]$mu -1.28*sd), 
                          ymax=plogis(val + repList[[bm]]$mu +1.28*sd)), alpha=0.3)
        scale_color_grey()
      }
      if(j==2){
        g[[j]][[i]] <- g[[j]][[i]] +
          geom_line(aes(x=x, y=val))+#, colour="black")) +
          geom_ribbon(aes(ymin=val-1.28*sd, ymax=val+1.28*sd), alpha=0.3)
        scale_color_grey()
      }
      
      if(myCat[j] == "re"){
        g[[j]][[i]] <- g[[j]][[i]] + ylab(" Marginal effect")
      }
      if(myCat[j] == "est"){
        g[[j]][[i]] <- g[[j]][[i]] + ylab("Surival")
        g[[j]][[i]] <- g[[j]][[i]] + geom_point(aes(x=x,y=obs))
      }
      g[[j]][[i]] <- g[[j]][[i]] + xlab(xlab[i])
      g[[j]][[i]] <- g[[j]][[i]] + theme(legend.position = "none")
      # print(g[[j]][[i]])
    }
  }
  
  # png("ggplot_marginal_effects_bestfit.png", height=640, width=480, res=100)
  
  ggarrange(g[[1]][[1]], g[[2]][[1]], 
            g[[1]][[2]], g[[2]][[2]], 
            g[[1]][[3]], g[[2]][[3]],
            labels = LETTERS[1:6],
            ncol = 2, nrow = 3)
  
  # dev.off()
  
}

