f_management_aggregate_comparison <- function(fit = fit,
                                              save_to_file = FALSE,
                                              width = 800,
                                              height = 400,
                                              res = 100,
                                              actions_to_include = NA,
                                              point_size = 20){
  library(ggplot2)
  
  proj <- t(as.list(fit$opt$SD, 'Est', report=TRUE)['proj'][[1]])
  projSD <- t(as.list(fit$opt$SD, 'Std', report=TRUE)['proj'][[1]])
  
  id_order <- (1:nrow(proj))[order(-proj[,1])]
  proj <- proj[id_order,]
  projSD <- projSD[id_order,]
  
  actionNames <- paste0 (fit$proj$grid[,'l'], ' mm',"\n",fit$proj$grid[,'j'], ' days')[]
  print(actionNames)
  
  df <- data.frame(proj = c(proj), sd = c(projSD))
  df$id <- rep(1:nrow(proj),ncol(proj))
  df$a <- rep(c("LGR","by-pass","aggregate"),each = nrow(proj))
  df$action <- rep(id_order,ncol(proj))
  
  df <- df[df$action %in% actions_to_include,]
  
  g <- ggplot(df[df$a!="aggregate",], aes(y=proj, x=as.factor(id), group=a)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Smolt to adult survival") +
    xlab("Management action") +
    geom_errorbar(aes(ymin=(proj-1.96*sd), ymax=(proj+1.96*sd), colour = a), width=.2,
                  position=position_dodge(0.5)) +
    ylim(0.006,0.0135) +
    geom_hline(yintercept = sum(fit$df$ns[fit$df$a==0])/sum(fit$df$nt[fit$df$a==0]), colour="#619CFF", size=1, alpha=0.5) +
    geom_hline(yintercept = sum(fit$df$ns[fit$df$a==1])/sum(fit$df$nt[fit$df$a==1]), colour="#F8766D", size=1, alpha=0.5) +
    geom_hline(yintercept = sum(fit$df$ns[])/sum(fit$df$nt[]), colour = "grey", size=3, alpha=0.5) +
    geom_point(aes(x=as.factor(id), y=(proj),colour=a),
               size=4,
               shape=15,
               position=position_dodge(0.5)) +
    theme(axis.text.x = element_text(angle = 0)) +
    scale_x_discrete(labels= actionNames[actions_to_include]) +
    guides(colour=guide_legend(title="Migration\npathway"))


  if(save_to_file){
    png("f_ggplot_management_aggregate_comparison.png",
        width=width, height = height,
        res=res,
        pointsize = point_size)
  }

  print(g)

  if(save_to_file) dev.off()
  return(g)
}

# dev.off()
