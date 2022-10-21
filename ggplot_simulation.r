

load("simOut.rData")

Action <- c("-7 days, +4 mm", "-7 days, -4 mm", "+7 days, -4 mm", "+7 days, +4 mm", "+7 days, 0 mm", "-7 days, 0 mm",  "0 days, +4 mm", "0 days, -4 mm", "0 days, 0 mm")

# tiff("ggplot_simulation.tiff", width=500, height=500)

par(mfrow=c(3,3))

# tmp <- aggregate(list(s=data$surv, t=data$total), by=list(data$t_i),sum)
# tmp <- tmp$s/tmp$t

for(i in 1:9){
  plot(c(t(na.omit(simOut[,,i,1]))),c(t(na.omit(simOut[,,i,2]))),
       pch=16,
       col="grey",
       cex=1.2,
       main= Action[i],
       ylab="With random effects",
       xlab="Wihtout random effects")
  segments(0,0,20,20)  
}

# dev.off()
