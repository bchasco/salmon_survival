f_ggplot_marginal_effects2 <- function(fit = fit,
                                      vars = 'mar_l',
                                      lab_hjust = 1,
                                      alpha_data = 0.1,
                                      include_data = FALSE,
                                      xNames = c('mar_l' = "Length (mm)", 'mar_j' = 'Calendar day', 'mar_y' = "Smolt year"),
                                      yRange = c('mar_l' = 0.08, 'mar_j' = 0.08, 'mar_y' = 0.06),
                                      qq = 1.64,
                                      height = 600,
                                      width = 400,
                                      point_size = 24,
                                      file = "f_ggplot_marginal_effects2.png",
                                      save_file = TRUE){

  library(ggplot2)
  library(ggpubr)
  

  plotlist <- list()
  icnt <- 1
  for(i in vars){
    sdEst <- as.data.frame(as.list(fit$opt$SD, "Est", report=TRUE)[i])
    ni <- nrow(sdEst)
    if(i!='mar_y'){
      sdEst <- sdEst %>%
        mutate(var = i, x = 1:ni + 84) %>%
        rename("PIT tagged" = !!names(.[1]),
               "bypass" = !!names(.[2])) %>%
        pivot_longer(cols = c("PIT tagged", "bypass")) %>%
        mutate(sd = c(t(as.data.frame(as.list(fit$opt$SD, "Std", report=TRUE)[i])))) 
    }else{
      sdEst <- sdEst %>%
        mutate(var = i, x = 1:ni + 1997) %>%
        rename("PIT tagged" = !!names(.[1]),
               "bypass" = !!names(.[2])) %>%
        pivot_longer(cols = c("PIT tagged", "bypass")) %>%
        mutate(sd = c(t(as.data.frame(as.list(fit$opt$SD, "Std", report=TRUE)[i])))) 
    }
    
    
    if(i == 'mar_j'){
      g <- ggplot(sdEst, aes(x=as.Date(x, origin = as.Date("2018-01-01")), 
                             y=exp(value), color = name)) +
        geom_line() + 
        geom_ribbon(aes(ymin = exp(value - 1.28 * sd), 
                        ymax = exp(value + 1.28 * sd),
                        fill = name),
                    alpha =0.5,
                    colour = NA) +
      ylab("") +
      # geom_text(aes(label = LETTERS[icnt],  x = as.Date(Inf, origin = as.Date("2018-01-01")), y = Inf), inherit.aes = FALSE, hjust = 0, vjust = 1) +
      guides(fill=guide_legend("LGR migration\npathway"), color = guide_legend("LGR migration\npathway")) +
      xlab(xNames[i])
    }
    
    if(i != 'mar_j'){
        g <- ggplot(sdEst, aes(x=x, y=exp(value), color = name)) +
        geom_line() + 
        geom_ribbon(aes(ymin = exp(value - 1.28 * sd), 
                        ymax = exp(value + 1.28 * sd),
                        fill = name),
                    alpha =0.5,
                    colour = NA) +
      ylab("") +
      # geom_text(aes(label = LETTERS[icnt],  x = -Inf, y = Inf), inherit.aes = FALSE, hjust = -1, vjust = 1) +
      guides(fill=guide_legend("LGR migration\npathway"), color = guide_legend("LGR migration\npathway")) +
      xlab(xNames[i])
    }
    
    if(i == 'mar_y'){
      g <- ggplot(sdEst, aes(x=x, y=exp(value), color = name)) +
        geom_line() + 
        geom_ribbon(aes(ymin = exp(value - 1.28 * sd), 
                        ymax = exp(value + 1.28 * sd),
                        fill = name),
                    alpha =0.5,
                    colour = NA) +
        ylab("") +
        guides(fill=guide_legend("LGR migration\npathway"), color = guide_legend("LGR migration\npathway")) +
        xlab(xNames[i])
      
    }
    
    plotlist[[icnt]] <- g + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    
    icnt <- icnt + 1
  }
  
  gp <- ggarrange(plotlist = plotlist,
                  hjust = lab_hjust,
                  labels = LETTERS[1:length(plotlist)],
                  ncol = 1, nrow = length(plotlist),
                  common.legend = TRUE,
                  legend = 'right')
  
  require(grid)
  p = annotate_figure(p = gp, left = textGrob("Survival", rot = 90, vjust = 1.5, gp = gpar(cex = 1.3)))
  
  if(save_file) png(file, height = height, width = width, pointsize = point_size)
  print(gp)
  if(save_file) dev.off()
  
  print(p)
  
  return(list(gp = gp,
         plotlist = plotlist,
         p = p))
}

# f_ggplot_marginal_effects(fit = fit,
#                           var = 'mar_l')

gp <- f_ggplot_marginal_effects2(fit = fit, vars = c('mar_l','mar_j','mar_y'), lab_hjust = 1)
