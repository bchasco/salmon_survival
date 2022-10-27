library(TMB)
library(INLA)
library(dplyr)
library(tidyr)

rm(list=ls())

source("f_raw_data.r") #Collect the raw data
source("f_mesh.r") #Build the mesh
source("f_map.r") #Build the mesh
source("f_TMB_data_pars.r") #Build the data and parameters lists
source("f_model.r") #Sned everything to a TMb object and estimate
source("f_cpp.r") #Compile new model
source("f_ggplot_interaction.r") #plot epsilon and or survival

f_cpp("v11_6") #link the library

fit <- f_model(raw_file = "C:/NOAA/LARGE_data/chin_survival_1998_2019_with_tagged_above.csv",
               rangeL = c(minL = 85, maxL = 125),
               rangeJ = c(minJ = 85, maxJ = 160),
               AR_flags = c(t_AR = 1, j_AR = 1, l_AR = 1),
               bypass_flags = c(t_bypass = 0, j_bypass = 0, 
                                l_bypass = 0, jlt_bypass = 0, 
                                mu_bypass = 1),
               re_flags = c(t_flag = 1, j_flag = 1, 
                            l_flag = 1, jlt_flag = 1),
               H_flag = 1,
               version = "v11_6",
               n_knots = 250,
               random = c("y_re", "l_re", 
                          "j_re", "z_jlt"),
               compare_AIC = TRUE,
               getsd = TRUE)

f_ggplot_interaction(fit = fit, bypass_cond = 1)

# source("table_AIC_comp.r")
# #If you decide to change something dramatic about the model
# clear_AIC_comp <- FALSE
# if(clear_AIC_comp){
#   AIC_comp <- list()
#   save(file="AIC_comp.rData",AIC_comp)  
# }
