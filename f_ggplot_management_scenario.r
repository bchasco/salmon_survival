f_ggplot_management_scenarios <- function(fit = fit,
                                   alpha = 0.6,
                                   n_sim = 1,
                                   save_to_file = TRUE,
                                   file = "f_ggplot_management_scenarios.png",
                                   action = 1){
  source("f_tau_sim.r")
  source('f_range_sim.r')
  library(ggpubr)
  
  g_range <- f_range_sim(fit = fit, 
                         action = action, 
                         alpha = alpha, 
                         n_sim = n_sim)
  
  g_tau <- f_tau_sim(fit = fit, 
                     action = action, 
                     alpha = alpha, 
                     n_sim = n_sim)  
  
  plotlist <- list(g_range = g_range,
                   g_tau = g_tau)
  
  p <- ggarrange(plotlist = plotlist,
                 labels = LETTERS[1:length(plotlist)],
                 ncol = 1)##length(plotlist))
  
  # Annotate the figure by adding a common labels
  annotate_figure(p,
                  bottom = text_grob("", color = "black",
                                     hjust = 1, x = 0.65, size = 11))
  
  print(p)
  if(save_to_file) png(file, height=400, width=500, pointsize = 12)
  print(p)
  if(save_to_file) dev.off()
}