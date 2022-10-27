f_ggplot_marginal_effects <- function(fit = fit,
                                      var = 'mar_l'){

  library(ggplot2)
  library(ggpubr)
  
  sdEst <- as.list(fit$opt$SD, "Est", report=TRUE)
  sdSd <- as.list(fit$opt$SD, "Std", report=TRUE)
  mar_re <- c("mar_l","mar_j","mar_y")
  
  x <- fit$df[,]
  if(var=='mar_l'){
    x <- x %>% group_by(var = l, var2 = j)%>% summarise(ns = sum(ns), 
                                                        nt = sum(nt))  #get rid of diff column
  }
  if(var=='mar_j'){
    x <- x %>% group_by(var = j, var2 = l)%>% summarise(ns = sum(ns), 
                                                        nt = sum(nt))  #get rid of diff column
  }
  if(var=='mar_y'){
    x <- x %>% group_by(var = y)%>% summarise(ns = sum(ns), 
                                            nt = sum(nt))  #get rid of diff column
  }

  re <- data.frame(x = 1:length(unlist(sdEst[var])), 
            y = unlist(sdEst[var]), 
            sd = unlist(sdSd[var]),
            mu = unlist(sdEst$mu)) 
  g <- ggplot(re, aes(x=x, y=plogis(y))) +
    geom_point(data = x, 
               aes(x = var - min(var) + 1, 
                   y = ns/nt, 
                   colour="blue",
                   alpha = 0.1), inherit.aes = FALSE) + 
    theme_bw() + 
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
          geom_line() +
          geom_ribbon(aes(ymin=plogis(y - 1.96*sd),
                          ymax=plogis(y + 1.96*sd)), alpha=0.5, size=1) +
  ylim(0,0.04)
        # scale_color_grey()
  print(g)
}

f_ggplot_marginal_effects(fit = fit,
                          var = 'mar_y')
