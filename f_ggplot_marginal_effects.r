f_ggplot_marginal_effects <- function(fit = fit,
                                      vars = 'mar_l',
                                      alpha = alpha,
                                      xNames = c('mar_l' = "Length (mm)", 'mar_j' = 'Calendar day', 'mar_y' = "Smolt year"),
                                      yRange = c('mar_l' = 0.08, 'mar_j' = 0.08, 'mar_y' = 0.06),
                                      qq = 1.64,
                                      height = 600,
                                      width = 400,
                                      file = "f_ggplot_marginal_effects.png",
                                      save_file = TRUE){

  library(ggplot2)
  library(ggpubr)
  

  plotlist <- list()
  icnt <- 1
  for(i in vars){
    sdEst <- as.list(fit$opt$SD, "Est", report=TRUE)[i]
    sdSd <- as.list(fit$opt$SD, "Std", report=TRUE)[i]
    x <- fit$df[,]
    if(i =='mar_l'){
      x <- x %>% 
        filter(a == 0) %>%
        group_by(var = l, j = j)%>% 
        summarise(ns = sum(ns),
                  nt = sum(nt))  #get rid of diff column
      myNames <- c("early","average","late")
    }
    if(i =='mar_j'){
      x <- x %>% 
        filter(a == 0) %>%
        group_by(var = j, l = l)%>% 
        summarise(ns = sum(ns),
                  nt = sum(nt))  #get rid of diff column
      myNames <- c("large","average","small")
    }
    if(i == 'mar_y'){
      x <- x %>%       
        filter(a == 0) %>%
        group_by(var = y)%>% 
        summarise(ns = sum(ns), nt = sum(nt))  #get rid of diff column
      myNames <- c("large/ \nearly","average","small/ \nlate")
    }
    
    re <- data.frame(x = rep(1:nrow(sdEst[[1]]),ncol(sdEst[[1]])),
                     Level = rep(myNames, each = nrow(sdEst[[1]])),
                     y = c(sdEst[[1]]), 
                     sd = c(sdSd[[1]])) 
    
    g <- ggplot(re, aes(x=x, y=plogis(y), group = Level)) +
      geom_point(data = x, 
                 aes(x = var - min(var) + 1, 
                     y = ns/nt, 
                     alpha = 0.1), 
                 show.legend = FALSE,
                 inherit.aes = FALSE) + 
      theme_bw() + 
      theme(panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            axis.line = element_line(colour = "black")) +
      geom_ribbon(aes(ymin=plogis(y - qq*sd),
                      ymax=plogis(y + qq*sd),
                      fill = Level),
                  alpha = 0.2) +
      geom_line(aes(colour = Level), size = 1) +
      guides(fill=guide_legend(title=""), colour=guide_legend(title="")) +
      xlab(xNames[i]) +
      ylim(0,yRange[i]) + 
      ylab("")
    plotlist[[icnt]] <- g
    icnt <- icnt + 1
  }
  gp <- ggarrange(plotlist = plotlist,
            labels = LETTERS[1:length(plotlist)],
            ncol = 1, nrow = length(plotlist))
  
  require(grid)
  annotate_figure(p = gp, left = textGrob("Survival", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
  if(save_file) png(file, height = height, width = width, pointsize = 16)
  print(gp)
  if(save_file) dev.off()
  return(gp)
}

# f_ggplot_marginal_effects(fit = fit,
#                           var = 'mar_l')
