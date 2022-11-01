#plot sim fixed effects
f_plot_sim_fixed_effects <- function(fits = list(Orig_data = fit, N_i_2X = fit_sim_size_2)){
  ex_names <- names(fits)
  varNames <- names(fits[[1]]$obj$par)
  varNames[varNames=="ln_H_input_jl"] <- c("lnH_1","lnH_2")
  
  df <- expand.grid(ex=ex_names,
                    var=varNames,
                    ex_i=1:nrow(fits[[1]]$bias_sim$sim_out), 
                    val = NA)

  mle <-   expand.grid(ex=ex_names,var=varNames,
                             ex_i=1:1, 
                             val = NA)
  
  for(i in ex_names){
    df$val[df$ex==i] <- c(t(fits[[i]]$bias_sim$sim_out))
    mle$val[mle$ex==i] <- c(fits[['Orig_data']]$opt$par)
  }
  print(ex_names)
  print(varNames)
  print(length(fits))
  
  library(ggplot2)
  g <- ggplot(df, aes(y = val, x = ex, colour = ex, fill = ex)) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    geom_boxplot(outlier.shape=8,
                 outlier.size=4, alpha = 0.5)+
    # geom_point(aes(colour = ex)) +
    geom_point(data = mle, aes(y = val, x = ex,size = 3), 
               shape=16, 
               colour = "black",  
               inherit.aes = FALSE,
               show.legend = FALSE) +
    facet_wrap(~var, scales = 'free') +
    
    
    
  
  print(g)

  return(df)
}

tmp <-f_plot_sim_fixed_effects()
