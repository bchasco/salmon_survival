load("fit_SD.rdata")
f_management_aggregate_comparison(fit = fit,
                                  save_to_file = TRUE,
                                  actions_to_include = c(1:9),
                                  width = 700, 
                                  height = 500)

f_management_annual_comparison(fit = fit,
                               save_to_file = TRUE,
                               actions_to_include = c(1:9),
                               width = 700,
                               height = 500)

f_ggplot_interaction(fit = fit,
                     bypass_cond = 1,
                     years_to_plot = c(1998:1999,2018:2019),
                     plot_type = c("both"),
                     gg_ncol = 1,
                     observed_max = 0.2,
                     include_data = TRUE,
                     save_to_file = TRUE)
# 
f_ggplot_marginal_effects(fit = fit,
                          vars = c('mar_j', 'mar_l', 'mar_y'),
                          qq = 01.28,
                          include_data = FALSE,
                          alpha_data = 0.05,
                          point_size = 50,
                          height = 500,
                          width = 300,
                          yRange = c('mar_l' = 0.045, 'mar_j' = 0.02, 'mar_y' = 0.035),
                          save_file = TRUE)

load("fit_noSD.rData")
tau_sim <- f_tau_sim(fit = fit, n_sim = 100, 
                     save_to_file = FALSE, 
                     res=100, 
                     seed = 200,
                     point_size = 12,
                     width = 400, 
                     height = 500)

range_sim <- f_range_sim(fit = fit, n_sim = 100, 
                         save_to_file = FALSE, 
                         res = 100, 
                         seed = 200,
                         point_size = 12,
                         width = 400, 
                         height = 500)
p <- ggarrange(plotlist = list(tau_sim,range_sim),
               labels = LETTERS[1:2],
               ncol = 2)
png(file = "Interaction_simulation.png", width = 500, height = 300)
print(p)
dev.off()
