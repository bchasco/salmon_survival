f_management_annual_comparison <- function(fit = fit,
                                              save_to_file = FALSE,
                                              actions_to_include = c(3)){
  library(ggplot2)

  proj <- (as.list(fit$opt$SD, 'Est', report=TRUE)['proj_y'][[1]])
  projSD <- (as.list(fit$opt$SD, 'Std', report=TRUE)['proj_y'][[1]])
  
  actionNames <- paste0 ('length ',fit$proj$grid[,'l'],
                                   "\n",'arrival ',
                                   fit$proj$grid[,'j'])
  
  df <- data.frame(proj = c(proj), sd = c(projSD))
  df$y <- rep(1:dim(proj)[3]+1997,each = dim(proj)[1] * length(actionNames))
  df$a <- rep(rep(c("LGR","by-pass","aggregate"),times = length(actionNames)), times = dim(proj)[3])
  df$action <- rep(rep(actionNames,each = dim(proj)[1]),dim(proj)[3])
  
  df <- df[df$action %in% actionNames[actions_to_include],]
  # for(m in unique(df$action)){
    # for(m in unique(df$action)[c(1:4,6:9,5)]){
    #   df$proj[df$action==m] <- 
    #     (df$proj[df$action==m] - df$proj[df$action==actionNames[5]])/df$proj[df$action==actionNames[5]]
    # }
  # }

  g <- ggplot(df[df$a!="aggregate",], aes(y=proj, x=y, group=a)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Percent change in SAR survival") +
    xlab("Migration year") +
    geom_ribbon(aes(ymin=(proj-1.96*sd),
                    ymax=(proj+1.96*sd),
                    fill = a),
                alpha = 0.3, colour = NA, show.legend = FALSE) +
    geom_hline(yintercept = 0) + 
    geom_line(aes(x=y, y=(proj),colour=a),
               size=0.5,
               shape=15) +
    theme(axis.text.x = element_text(angle = 0)) +
    facet_wrap(~action, scales="free") + 
    # scale_x_discrete(labels= actionNames[id_order]) + 
    guides(colour=guide_legend(title="Migration\npathway"))
  
  
  if(save_to_file){
    png("f_ggplot_management_annual_comparison.png",
        width=500, height = 500,
        res=100)
  }

  print(g)

if(save_to_file) dev.off()
}
