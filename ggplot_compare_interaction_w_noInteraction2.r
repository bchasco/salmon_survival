mu <- apply(simOut,c(2,3,4),mean)
sd <- apply(simOut,c(2,3,4),sd)

sim <- data.frame(x = c(simOut[,,,2]), 
                  nox = c(simOut[,,,1]),
                  m = rep(1:9, each=19*dim(simOut)[1]))
library(ggplot2)

p <- ggplot(data=sim, aes(x=x,y=nox)) + 
  geom_point() + 
  facet_wrap(vars(m), nrow=3) +
  geom_abline(slope=1,intercept = 0)
print(p)
