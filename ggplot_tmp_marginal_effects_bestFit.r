library(ggplot2)
library(ggpubr)


sdEst <- sd$value[names(sd$value)=="l_re"]#sdList$est[[bm]]
sdSd <- sd$sd[names(sd$value)=="l_re"]#sdList$sd[[bm]]
ai <- rep(c("Detector","By-pass"), each=length(sdSd)/2)
mu <- rep(c(rep$mu, rep$mu ), each=length(sdSd)/2)
# x <- rep(1998:2019,2)
# x <- rep(minJ:maxJ,2)
x <- rep(minL:maxL,2)

df <- data.frame(mu = sdEst,sdSd,ai,x = x)
g <- ggplot(df,aes(x = x,y = (mu),color=as.factor(ai))) +
  geom_ribbon(aes(ymin = (mu - sdSd*1.28),
                  ymax = (mu + sdSd*1.28),
                  fill = as.factor(ai)),
              alpha=0.2, colour = NA) +
  geom_line(aes(color=as.factor(ai)), size = 1) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  
print(g)

