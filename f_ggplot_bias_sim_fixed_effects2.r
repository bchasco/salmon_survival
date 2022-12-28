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
  mle$val <- fit$opt$par
  
  library(ggplot2)
  g <- ggplot(df, aes(y = val)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_boxplot(outlier.shape=8,
                 outlier.size=4,
                 fill = "grey",
                 alpha = 0.5)+
    facet_wrap(~var, scales = 'free') +
    geom_point(data = mle, aes(y = val, x = 0, size = 3),
               shape=16,
               colour = "red",
               alpha = 0.5,
               inherit.aes = FALSE,
               show.legend = FALSE) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
    
    
    
    
  
  print(g)

  return(df)
}

tmp <- f_ggplot_bias_sim_fixed_effects2(fit = fit)
