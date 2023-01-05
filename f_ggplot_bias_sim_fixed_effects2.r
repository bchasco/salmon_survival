#plot sim fixed effects
f_ggplot_bias_sim_fixed_effects2 <- function(fit = fit){
  
  varNames <- names(fit$obj$par)
  if('ln_H_input_jl'%in%varNames)
    varNames[varNames=="ln_H_input_jl"] <- c("lnH_1","lnH_2")
  
  df <- expand.grid(var=varNames,
                    sim_i=1:nrow(fit$bias_sim$sim_out))
  df$val <- c(t(fit$bias_sim$sim_out))

  mle <- expand.grid(var=varNames,
                       sim_i=1:1)
  mle$val <-fit$opt$SD$par.fixed
  # mle$val[1] <- -4.5
  mle$SD <- sqrt(diag(fit$opt$SD$cov.fixed))
  
  pd <- position_dodge(0.1) # move them .05 to the left and right
  
  library(ggplot2)
  g <- ggplot(df, aes(y = val)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_boxplot(outlier.shape=NA,
                 outlier.size=NA,
                 fill = "grey",
                 outlier.colour=NA,
                 alpha = 0.5)+
    facet_wrap(~var, scales = 'free') +
    geom_point(data = mle, aes(y = val, x = 0.1, size = 3),
               shape=16,
               colour = "red",
               alpha = 0.5,
               inherit.aes = FALSE,
               show.legend = FALSE) +
    geom_errorbar(data=mle,
                  aes(y = val, x = 0.1, ymin=val-1.28*SD, ymax=val+1.28*SD), 
                  # inherit.aes = FALSE,
                  colour = 'red',
                  width=.1, 
                  position=pd) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
    
    
    
  ylim1 = boxplot.stats(df$val)$stats[c(1, 5)]
  # g = g + coord_cartesian(ylim = ylim1*1.05)
  print(g)
  print(ylim1)
  return(df)
}

# load("bias_sim_re_y_j_l_AR1_j_l_bypassmu.rData")
# load("bias_sim_re_j_AR1_j_bypass_none.rData")
# load("bias_sim_re_y_AR1_y_bypass_none.rData")
# load("bias_sim_re_l_iid_l_bypass_none.rData")
# load("bias_sim_re_y_j_l_AR1_j_l_bypassmu.rData")
# load("bias_sim_re_y_RW_y_bypass_none.rData")
# load("bias_sim_Model_55.rData")
# tmp <- f_ggplot_bias_sim_fixed_effects2(fit = fit)
