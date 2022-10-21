# source("anisotropy_plot.r")
library(ggpubr)

library(RandomFields)
library(raster)
library(RANN)
library(ggplot2)

set.seed(1)
load("repList.rData")
bm <- 2
myList <- obj$simulate(complete=TRUE)
orgRep <- repList[[bm]]

nt <- length(parameters$y_re)
DF <- data.frame(do.call(rbind,replicate(nt, vizloc_xy, simplify=FALSE)),
                 rep(1:nt+1997,each=nrow(vizloc_xy)))

DF$re <- c(apply(myList$ln_zexp_sp_jlt,2,function(x){return(x[knots_xy$nn.idx])}))
yr_ef <- data.frame(y=1998:2016,ef=plogis(myList$mu+myList$y_re*myList$sig_y)) #This mu is incorrect
j_ef <- data.frame(j=min(df2$j):max(df2$j),ef=myList$j_re*myList$sig_j) #This mu is incorrect
l_ef <- data.frame(l=minL:maxL,ef=myList$l_re*myList$sig_l) #This mu is incorrect

tiff("ggplot_y_re_simulation.tiff", height=400, width=400)

plot(yr_ef, type="l",
     las=1,
     ylab="Survival",
     xlab="Migration year")
points(data.frame(y=1998:2016,ef=plogis(orgRep$mu+orgRep$y_re*orgRep$sig_y)),
       pch=16, col="grey", cex=1.5) #This mu is incorrect

dev.off()
