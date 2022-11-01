# rm(list=ls())

library(TMB)
library(INLA)
library(dplyr)
library(tidyr)
library(DHARMa)
library(raster)
library(RANN)


source("f_cpp.r") #Compile new model
source("f_raw_data.r") #Collect the raw data
source("f_management_projection.r") #Compile new model
source("f_mesh.r") #Build the mesh
source("f_map.r") #Build the mesh
source("f_TMB_data_pars.r") #Build the data and parameters lists
source("f_model.r") #Sned everything to a TMb object and estimate
source("f_bias_sim.r")

source("f_ggplot_interaction.r") #plot epsilon and or survival
source("f_ggplot_marginal_effects.r") #plot epsilon and or survival
source("f_management_aggregate_comparison.r")
source("f_management_annual_comparison.r")

f_cpp("v11_6") #link the library

file <- "C:/NOAA/LARGE_data/chin_survival_1998_2019_with_tagged_above.csv"

sys <- Sys.time()

fit <- f_model(raw_file = file,
               
               #Data trimming
               rangeL = c(minL = 85, maxL = 125),
               rangeJ = c(minJ = 85, maxJ = 160),
               
               #Number of knots
               n_knots = 250,
               
               #Management actions for projecting
               #j in days, l in mm
               m_proj = data.frame('j' = c(-7,0,7),
                                   'l' = c(-4,0,4)),
               
               #Random effects (0 = leave out, 1 = include)
               re_flags = c(t_flag = 0, j_flag = 0, 
                            l_flag = 0, jlt_flag = 0),
               
               #R.E. statistical model (0 = i.i.d, 1 = AR1, 2 = R.W.) 
               AR_flags = c(t_AR = 1, j_AR = 1, l_AR = 1),
               
               #Bypass offset (0 = FALSE, 1 = TRUE)
               bypass_flags = c(t_bypass = 0, j_bypass = 0, 
                                l_bypass = 0, jlt_bypass = 0, 
                                mu_bypass = 1),
               
               #TMB (random or fixed)
               random = c("y_re", "l_re", 
                          "j_re", "z_jlt"),
               
               H_flag = 1,
               version = "v11_6",
               compare_AIC = FALSE,
               getsd = FALSE,
               save_to_file = FALSE,
               save_file = "fit.rData",
               DHARMa_sim = FALSE,
               bias_sim = TRUE,
               bias_sim_n = 10,
               sim_size = 1,
               proj_sim = FALSE,
               proj_H = c(0.,0.,))

# f_management_aggregate_comparison(fit = fit,
#                                   save_to_file = FALSE)
# 
# f_management_annual_comparison(fit = fit,
#                                save_to_file = FALSE,
#                                actions_to_include = c(1,3,7,9))
# 
# f_ggplot_interaction(fit = fit,
#                      bypass_cond = 1,
#                      years_to_plot = 1998:2019,
#                      plot_type = "Survival")
# 
# f_ggplot_marginal_effects(fit = fit,
#                           vars = c('mar_j'),
#                           qq = 1.96,
#                           height = 400,
#                           width = 400,
#                           yRange = c('mar_l' = 0.08, 'mar_j' = 0.08, 'mar_y' = 0.08),
#                           save_file = TRUE,
#                           alpha = alpha)

# source("table_AIC_comp.r")
# #If you decide to change something dramatic about the model
# clear_AIC_comp <- FALSE
# if(clear_AIC_comp){
#   AIC_comp <- list()
#   save(file="AIC_comp.rData",AIC_comp)  
# }

print(Sys.time()-sys)