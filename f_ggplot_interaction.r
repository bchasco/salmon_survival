f_ggplot_interaction <- function(fit = fit,
                                 years_to_plot = c(1999,2019),
                                 file = "f_ggplot_interaction.png",
                                 save_to_file = FALSE,
                                 gg_ncol=4,
                                 arrange_ncol = 2,
                                 bypass_cond = 1,
                                 plot_type = "Survival"){

  library(ggpubr)
  library(RandomFields)
  library(raster)
  library(RANN)
  library(ggplot2)
  

  rep <- fit$obj$rep
  df <- fit$df
  mesh <- fit$mesh$mesh_jl
  data <- fit$obj$env$data
  
  # Generate grid to visualize density
  x_seq <- seq(min(df$j),max(df$j),by=1) #Sequence over days
  y_seq <- seq(min(df$l),max(df$l),by=1) #Sequence over lengths
  vizloc_xy = expand.grid( x=x_seq, y=y_seq ) #Field of length and day
  
  knots_xy = nn2(data=mesh$loc[,1:2], query=vizloc_xy, k=1) #Map the data to the knots using nearest-neighbor

  nt <- length(rep$y_re)
  DF <- data.frame(do.call(rbind,replicate(nt, vizloc_xy, simplify=FALSE)),
                   rep(1:nt+1997,each=nrow(vizloc_xy)))
  
  DF$re <- c(apply(rep$z_jlt[,bypass_cond,],2,function(x){return(x[knots_xy$nn.idx])}))
  yr_ef <- data.frame(y=1998:2019,
                      ef=rep$y_re*rep$sig_y) #This mu is incorrect
  j_ef <- data.frame(j=min(df$j):max(df$j),
                     ef=rep$j_re*rep$sig_j) #This mu is incorrect
  l_ef <- data.frame(l=min(df$l):max(df$l),
                     ef=rep$l_re*rep$sig_l) #This mu is incorrect
  
  names(DF) <- c("x","y","tt",'R.E.')
  DF$mu <- unlist(rep['mu'])
  DF$yr <- yr_ef$ef[match(DF$tt,yr_ef$y)]
  DF$j <- j_ef$ef[match(round(DF$x),j_ef$j)]
  DF$l <- l_ef$ef[match(round(DF$y),l_ef$l)]
  DF$Survival <- plogis(DF$mu+DF$yr+DF$j+DF$l+DF$R.E.)
  
  # DF2 <- data.frame(x=data$j_i+min(df$j), 
  #                   y=data$l_i+min(df$l), 
  #                   tt=data$t_i+1998, 
  #                   Observed=data$surv/data$total)
  
  sampleYears <- years_to_plot
  breaks <- c(-1,0,1)
  p1 <- ggplot(data = DF[DF$tt%in%sampleYears,], 
               aes(x = as.Date(x, origin = as.Date("2018-01-01")), 
                   y = y,
                   fill = R.E.)) +
    geom_raster()+
    facet_wrap(~tt, ncol=gg_ncol) +
    ylab('Length (mm)') +
    xlab('') +
    scale_fill_gradientn(colours=c("grey","blue", "red")) +
    # scale_x_continuous(expand = c(0,0)) + 
    scale_x_date(limits=range, expand =c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position="right")
  p2 <- ggplot(data = DF[DF$tt%in%sampleYears,], aes(x = as.Date(x, origin = as.Date("2018-01-01")), y = y, 
                                                     fill = Survival)) +
    geom_raster()+
    facet_wrap(~tt, ncol=gg_ncol) +
    ylab('Length (mm)') +
    xlab('') +
    scale_fill_gradientn(colours=c("grey","blue", "red")) +
    # scale_x_continuous(expand = c(0,0)) + 
    scale_x_date(limits=range, expand =c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position="right")
  
  # DF2$Observed[DF2$Observed>0.2] <- 0.2
  # p3 <- ggplot(data = DF2[DF2$tt%in%sampleYears,], aes(x = as.Date(x, origin = as.Date("2018-01-01")), y = y, 
  #                                                      fill = Observed)) +
  #   geom_raster()+
  #   facet_wrap(~tt, ncol=1) +
  #   ylab('Length (mm)') +
  #   xlab('') +
  #   scale_fill_gradientn(colours=c("grey","blue", "red")) +
  #   # scale_x_continuous(expand = c(0,0)) + 
  #   theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
  #                      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  #   scale_x_date(limits=range, expand =c(0,0)) +
  #   scale_y_continuous(expand = c(0,0)) +
  #   theme(legend.position="right")
  
  p <- ggarrange(p1, p2,
                 labels = LETTERS[1:2],
                 ncol = arrange_ncol)
  
  # Annotate the figure by adding a common labels
  annotate_figure(p,
                  bottom = text_grob("Migration date past Lower Granite Dam", color = "black",
                                     hjust = 1, x = 0.65, size = 11))

  if(save_to_file) tiff(file, height=820, width=820)
  if(plot_type == "RE") print(p1)
  if(plot_type == "Survival") print(p2)
  if(plot_type == "both") print(p)
  if(save_to_file) dev.off()
  
}