library(tidyr)
library(ggpubr)
library(RandomFields)
library(raster)
library(RANN)
library(ggplot2)
# Generate grid to visualize density
vizloc_xy <- expand.grid(j_i = seq(min(range(data$j_i)),
                      max(range(data$j_i)),
                      1)+minJ,
            l_i = seq(min(range(data$l_i)),
                      max(range(data$l_i)),
                      1)+minL)

h_seq <- rep(apply(hh,1,function(x){return(paste(round(x,2),collapse = ", "))}), 
             each = length(jj)*length(kk)*nrow(vizloc_xy))
j_seq <- rep(rep(jj, 
                 each = nrow(hh)*nrow(vizloc_xy)),
             times=length(kk))
k_seq <- rep(kk, 
             each = nrow(vizloc_xy), 
             times =nrow(hh)*length(jj))

n_per <- nrow(hh)*length(jj)*length(kk)

DF <- data.frame(do.call(rbind,
                         replicate(n_per, 
                                   vizloc_xy, 
                                   simplify=FALSE)),
                 rep(1:n_per,each=nrow(vizloc_xy)))

knots_xy = nn2(data=loc_xy_jl, query=vizloc_xy, k=1) #Map the data to the knots using nearest-neighbor

DF$h <- paste(h_seq,'H')
DF$j <- paste(j_seq,'tau')
DF$k <- paste(round(sqrt(8.0) / exp(kk),1),'range, ',round(exp(kk),2),'kappa')

RE <- NA
for(i in 1:n_per){
  RE <- c(RE,
          simRep_2X2$sim_i[[1]][[i]]$ln_zexp_sp_jlt[knots_xy$nn.idx,1])
}
RE <- na.omit(RE)  
DF$re <- RE

names(DF) <- c('x','y','s_i','H','tau','kappa','R.E.')

p1 <- ggplot(data = DF[,], aes(x = as.Date(x, origin = as.Date("2018-01-01")), 
                               y = y, 
                               fill = R.E.)) +
  facet_wrap(~H)+
  geom_raster()+
  ylab('Length (mm)') +
  xlab('Day') +
  scale_fill_gradientn(colours=c("grey","blue", "red")) 
print(p1)
