# source("anisotropy_plot.r")

library(ggpubr)
library(RandomFields)
library(raster)
library(RANN)
library(ggplot2)

# load("repList.rData")
# load("sdList.rData")
# bm <- 3
# myList <- repList[[bm]]

myList <- rep

# Generate grid to visualize density
x_seq <- seq(min(df2$j),max(df2$j),by=1) #Sequence over days
y_seq <- seq(min(df2$len),max(df2$len),by=1) #Sequence over lengths
vizloc_xy = expand.grid( x=x_seq, y=y_seq ) #Field of length and day

knots_xy = nn2(data=loc_xy_jl, query=vizloc_xy, k=1) #Map the data to the knots using nearest-neighbor
knots_obs = nn2(data=loc_xy_jl, query=cbind(data$j_i+minJ,data$l_i+minL), k=1) #Map the data to the knots using nearest-neighbor

nt <- length(parameters$y_re)
DF <- data.frame(do.call(rbind,replicate(nt, vizloc_xy, simplify=FALSE)),
                 rep(1:nt+1997,each=nrow(vizloc_xy)))

DF$re <- c(apply(myList$ln_zexp_sp_jlt,2,function(x){return(x[knots_xy$nn.idx])}))
yr_ef <- data.frame(y=1998:2019,ef=myList$y_re*myList$sig_y) #This mu is incorrect
j_ef <- data.frame(j=min(df2$j):max(df2$j),ef=myList$j_re*myList$sig_j) #This mu is incorrect
l_ef <- data.frame(l=minL:maxL,ef=myList$l_re*myList$sig_l) #This mu is incorrect

names(DF) <- c("x","y","tt",'R.E.')
DF$mu <- opt$par['mu']
DF$yr <- yr_ef$ef[match(DF$tt,yr_ef$y)]
DF$j <- j_ef$ef[match(round(DF$x),j_ef$j)]
DF$l <- l_ef$ef[match(round(DF$y),l_ef$l)]
DF$Survival <- plogis(DF$mu+DF$yr+DF$j+DF$l+DF$R.E.)

DF2 <- data.frame(x=data$j_i+minJ, y=data$l_i+minL, tt=data$t_i+1998, Observed=data$surv/data$total)

sampleYears <- c(1999:2019)#c(1998:2000,2017:2019)
breaks <- c(-1,0,1)
p1 <- ggplot(data = DF[DF$tt%in%sampleYears,], aes(x = as.Date(x, origin = as.Date("2018-01-01")), y = y, 
                               fill = R.E.)) +
                               geom_raster()+
  facet_wrap(~tt, ncol=round(sqrt(length(sampleYears)),0)) +
  ylab('Length (mm)') +
  xlab('') +
  scale_fill_gradientn(colours=c("grey","blue", "red")) +
  # scale_x_continuous(expand = c(0,0)) + 
  scale_x_date(limits=range, expand =c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme(legend.position="right")

print(p1)
