ggplot_interaction_function <- function(obj,
                                        sampleYears = 1998:2019,
                                        l_range = c(minL,maxL),
                                        j_range = c(minJ,maxJ),
                                        ncol = 4,
                                        grp = "Observed",
                                        simulated = TRUE){
  
  library(ggpubr)
  library(RandomFields)
  library(raster)
  library(RANN)
  library(ggplot2)
  
  rep <- sim_obj$report()

  # Generate grid to visualize density
  x_seq <- seq(j_range[1],
               j_range[2],
               by=1) #Sequence over days
  y_seq <- seq(l_range[1],
               l_range[2],
               by=1) #Sequence over lengths
  
  vizloc_xy = expand.grid(x = x_seq, 
                          y = y_seq ) #Field of length and day
  
  knots_xy = nn2(data = obj$mesh$loc[,1:2],
                 query = vizloc_xy,
                 k = 1) #Map the data to the knots using nearest-neighbor
  
  nt <- obj$env$.data$n_t
  
  #2D locations by year
  DF <- data.frame(do.call(rbind,
                           replicate(nt,
                                     vizloc_xy,
                                     simplify=FALSE)),
                   rep(1:nt+1997,
                       each=nrow(vizloc_xy))
                   )
  
  #2d random effects
  if(simulated){
    DF$re <- c(apply(sim_i$ln_zexp_sp_jlt,
                     2,
                     function(x){return(x[knots_xy$nn.idx])}))
  }else{
    DF$re <- c(apply(rep$ln_zexp_sp_jlt,
                     2,
                     function(x){return(x[knots_xy$nn.idx])}))
    
  }
  
  #rename things
  names(DF) <- c("x","y","tt",'R.E.')
  
  #AR1 indexes and random effects
  yr_ef <- data.frame(y = 1998:2019,
                      ef = rep$y_re) #This mu is incorrect
  j_ef <- data.frame(j = j_range[1]:j_range[2],
                     ef = rep$j_re) #This mu is incorrect
  l_ef <- data.frame(l = l_range[1]:l_range[2],
                     ef = rep$l_re) #This mu is incorrect
  
  #Add random effects to data.frame
  DF$mu <- opt$par['mu']
  DF$yr <- yr_ef$ef[match(DF$tt,yr_ef$y)]
  DF$j_ef <- j_ef$ef[match(round(DF$x),j_ef$j)]
  DF$l_ef <- l_ef$ef[match(round(DF$y),l_ef$l)]
  DF$Survival <- plogis(DF$mu+DF$yr+DF$j+DF$l+DF$R.E.)
  
  #Separate array for observations
  DF2 <- data.frame(x = obj$env$data$j_i+j_range[1], 
                    y = obj$env$data$l_i+l_range[1], 
                    tt = obj$env$data$t_i+1998, 
                    Observed = obj$env$data$surv/obj$env$data$total,
                    Total = obj$env$data$total)
  
  #Select the correct array
  ifelse(grp=="Observed" | grp=="Total",
         plotDF <- DF2,
         plotDF <- DF)
  
  #Create a raster plot of 2D results
  p1 <- ggplot(data = plotDF[plotDF$tt%in%sampleYears,],
               aes(x = as.Date(x, origin = as.Date("2018-01-01")),
                   y = y,
                   fill = plotDF[,grp]))+
    geom_raster() +
    facet_wrap(~tt, ncol = ncol) +
    ylab('Length (mm)') +
    xlab('') +
    scale_fill_gradientn(colours=c("grey","blue", "red")) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw() + 
    theme(panel.border = element_blank(), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), 
                       axis.line = element_line(colour = "black")) +
    scale_x_date(limits=range, expand =c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position="right") +
    labs(fill = grp)
    
    print(p1)
  return()
}

ggplot_interaction_function(sim_obj,
                            grp = "R.E.",
                            simulated=FALSE)
