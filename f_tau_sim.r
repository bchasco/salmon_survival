f_tau_sim <- function(fit = fit,
                        n_sim = 20,
                        seed = 100,
                        save_to_file = TRUE,
                      xlim = c(-100,150),
                      height = 300,
                      width = 300,
                      alpha = 0.6,
                      action = 1,
                      point_size = 16,
                      show.legend = TRUE,
                      res = 300){
  
  o.tau <- fit$obj$env$last.par['log_tau_jl2'] 
  sig <- 1/exp(o.tau)
  sig <- sig*c(0.5,1,2)
  proj_tau <- log(1/sig)
  
  set.seed(seed)
  proj_sim_n <- n_sim
  proj_sim <- array(NA, c(proj_sim_n,3,c(dim(fit$obj$rep$proj_y))))
  jcnt <- 1
  for(j in proj_tau){
    for(i in 1:proj_sim_n){
      fit$obj$env$data$proj_sim <- 1
      fit$obj$env$last.par[which(names(fit$obj$env$last.par) == 'log_tau_jl2')] <- j
      proj_sim[i,jcnt,,,] <- fit$obj$simulate(complete=TRUE)$proj_y
    }
    jcnt <- jcnt + 1
  }
  
  x0.5 <- round(sum(na.omit(c(proj_sim[,1,1:2,action,])<0))/
                  length(na.omit(c(proj_sim[,1,1:2,action,])))*100,1)
  df <- data.frame(val = c(proj_sim[,1,1:2,action,]),
                   var=rep(paste0("0.5x (",x0.5,"%)"),
                           length(c(proj_sim[,1,1:2,action,]))))
  
  x1 <- round(sum(na.omit(c(proj_sim[,2,1:2,action,])<0))/
                length(na.omit(c(proj_sim[,2,1:2,action,])))*100,1)
  df <- rbind(df,
              data.frame(val = c(proj_sim[,2,1:2,action,]),
                         var=rep(paste0("1x (",x1,"%)"),
                                 length(c(proj_sim[,2,1:2,action,])))))
  
  x2 <- round(sum(na.omit(c(proj_sim[,3,1:2,action,])<0))/
                length(na.omit(c(proj_sim[,3,1:2,action,])))*100,1)
  
  df <- rbind(df,
              data.frame(val = c(proj_sim[,3,1:2,action,]),
                         var=rep(paste0("2x (",x2,"%)"),
                                 length(c(proj_sim[,3,1:2,action,])))))
  
  g <- ggplot(df, aes(x=val)) +
    # facet_wrap(~var, ncol = 1, scales = "free") +
    geom_histogram(aes(fill = var), 
                   alpha=alpha, 
                   position="identity", 
                   show.legend = show.legend) +
    xlab("Percent change in annual survival") + 
    ylab("Frequency") +
    xlim(xlim) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(panel.background = element_rect(fill = NA)) + 
    geom_vline(xintercept = 0) + 
    guides(fill=guide_legend(title="Variance"))
  
  if(save_to_file){
    png(file = "f_tau_sim.png",
        width = width,
        height = height,
        pointsize = point_size,
        res = res)
    print(g)
    dev.off()
  }
  
  print(g)
  
  fit$obj$env$last.par['log_tau_jl2'] <- o.tau

  return(g)  
}

  
