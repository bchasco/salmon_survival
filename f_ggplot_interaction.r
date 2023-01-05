f_ggplot_interaction <- function(fit = fit,
                                 years_to_plot = c(1999,2019),
                                 file = "f_ggplot_interaction.png",
                                 save_to_file = FALSE,
                                 gg_ncol=2,
                                 arrange_ncol = 2,
                                 bypass_cond = 2,
                                 observed_max = 0.2,
                                 plot_type = "Survival",
                                 include_data = TRUE){

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

  nt <- nrow(rep$y_re)
  print(nt)
  print(1:nt+1997)
  DF <- data.frame(do.call(rbind,replicate(nt, vizloc_xy, simplify=FALSE)),
                   rep(1:nt+1997,each=nrow(vizloc_xy)))
  
  DF$re <- c(apply(rep$z_jlt[,1,],2,function(x){return(x[knots_xy$nn.idx])}))
  
  yr_ef <- data.frame(y=1998:2019,
                      ef=rep$y_re) #This mu is incorrect
  print(dim(DF))
  j_ef <- data.frame(j=min(df$j):max(df$j),
                     ef=rep$j_re) #This mu is incorrect
  l_ef <- data.frame(l=min(df$l):max(df$l),
                     ef=rep$l_re) #This mu is incorrect
  
  names(DF) <- c("x","y","tt",'R.E.')
  DF$mu <- unlist(rep['mu'])
  DF$yr <- yr_ef[match(DF$tt,yr_ef$y),2:ncol(yr_ef)]
  DF$j <- j_ef[match(round(DF$x),j_ef$j),2:ncol(j_ef)]
  DF$l <- l_ef[match(round(DF$y),l_ef$l),2:ncol(l_ef)]
  DF$noRE <- (DF$mu+as.matrix(DF$yr)+as.matrix(DF$j)+as.matrix(DF$l))
  DF$Survival <- (DF$mu+as.matrix(DF$yr)+as.matrix(DF$j)+as.matrix(DF$l) + rep(DF$R.E.))

  DF2 <- data.frame(x=data$j_i+min(df$j),
                    y=data$l_i+min(df$l),
                    a = data$a,
                    tt=data$t_i+1998,
                    surv = data$surv,
                    total = data$total)
  DF2 <- DF2 %>%
    group_by(x,y,tt) %>%  #get rid of diff column
    summarise(surv = sum(surv),
              total = sum(total)) %>% #survivors and sample size
    mutate(Observed = surv/total)

  sampleYears <- years_to_plot
  breaks <- c(-200,0,200)
  p1 <- ggplot(data = DF[DF$tt%in%sampleYears,],
               aes(x = as.Date(x, origin = as.Date("2018-01-01")),
                   y = y,
                   fill = plogis(Survival[,bypass_cond])/plogis(noRE[,bypass_cond])*100)) +
    geom_raster()+
    facet_wrap(~tt, ncol=gg_ncol) +
    ylab('Length (mm)') +
    xlab('') +
    scale_fill_gradientn(colours=c("blue","white", "red"),
                         guide = guide_legend(reverse = FALSE),
                         values = scales::rescale(c(0,100,300)),
                         breaks = c(0,100,300),
                         labels = c(0,100,300),
                         limits = c(0,300)) +
    # scale_x_continuous(expand = c(0,0)) +
    scale_x_date(limits=range, expand =c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position="right") +
    guides(fill=guide_legend(title="Interaction\neffect"))
  
  
  # print(head(DF[DF$tt%in%sampleYears,]))
  max_surv <- 0.0800#max(plogis(DF[DF$tt%in%sampleYears,]$Survival[,bypass_cond]))
  min_surv <- 0.0#min(plogis(DF[DF$tt%in%sampleYears,]$Survival[,bypass_cond]))
  mid_surv <- exp(fit$opt$SD$value[names(fit$opt$SD$value)=='tab_mu_a0'])
  breaks <- c(min_surv, mid_surv, max_surv)
  print(breaks)
  p2 <- ggplot(data = DF[DF$tt%in%sampleYears,], aes(x = as.Date(x, origin = as.Date("2018-01-01")), y = y,
                                                     fill = plogis(Survival[,bypass_cond]))) +
    geom_raster()+
    facet_wrap(~tt, ncol=gg_ncol) +
    ylab('Length (mm)') +
    xlab('') +
    scale_fill_gradientn(colours=c("blue","white","red"),
                         values = scales::rescale(breaks),
                         labels = round(breaks,4),
                         limits=c(0,max_surv),
                         breaks = breaks) +
    # scale_colour_gradient2(low = "blue", 
    #                        mid = "white", 
    #                        high = "red", 
    #                        midpoint = plogis(fit$obj$rep$mu + fit$obj$rep$bypass * (bypass_cond - 1)),
    #                        breaks = c(0.0,mid_surv,0.03)) +
    # scale_x_continuous(expand = c(0,0)) +
    scale_x_date(limits=range, expand =c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position="right") +
    guides(fill=guide_legend(title="Predicted\nSurvival"))


  DF2$Observed[DF2$Observed>=observed_max] <- observed_max
  mid_surv <- sum(DF2$surv)/sum(DF2$total)
  breaks <- c(0,mid_surv, observed_max)
  p3 <- ggplot(data = DF2[DF2$tt%in%sampleYears,], aes(x = as.Date(x, origin = as.Date("2018-01-01")), y = y,
                                                       fill = Observed)) +
    geom_raster()+
    facet_wrap(~tt, ncol=gg_ncol) +
    ylab('Length (mm)') +
    xlab('') +
    scale_fill_gradientn(colours=c("blue","white","red"),
                         values = scales::rescale(breaks),
                         labels = c("",round(breaks[2:3],4)),
                         limits=c(0,observed_max),
                         breaks = breaks,
                         guide = guide_legend(reverse = FALSE)) +
    # scale_x_continuous(expand = c(0,0)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(panel.background = element_rect(fill = 'lightgrey')) +
    scale_x_date(limits=range, expand =c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(legend.position="right")

  if(plot_type == "RE"){
    plotlist <- p1
    if(include_data)
      plotlist <- list(p1,p3)
  }

  if(plot_type == "Survival"){
    plotlist <- p2
    if(include_data)
      plotlist <- list(p2,p3)
  }

  if(plot_type == "both"){
    plotlist <- list(p1,p2)
    if(include_data)
      plotlist <- list(p1,p2,p3)
  }

  # p <- ggarrange(plotlist = plotlist,
  #                labels = LETTERS[1:length(plotlist)],
  #                ncol = length(plotlist))
  # 
  # # Annotate the figure by adding a common labels
  # annotate_figure(p,
  #                 bottom = text_grob("Migration date past Lower Granite Dam", color = "black",
  #                                    hjust = 1, x = 0.65, size = 11))
  # 
  print(p3)
  if(save_to_file) png(file, height=500, width=700, pointsize = 12)
  print(p3)
  if(save_to_file) dev.off()
  
  # print(p1)
  # print(p2)
  
  return(list(DF = DF
              ,DF2 = DF2
              ,yr_ef = yr_ef
              ,l_ef = l_ef
              ,j_ef = j_ef
              # ,DF2 = DF2
              ))
  
}
