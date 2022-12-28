rm(list=ls())

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
source("f_bias_sim.r") #Function for loop over bias estimates
# source("f_proj_sim.r") #Function for loop over bias estimates

source("f_ggplot_interaction.r") #plot epsilon and or survival
source("f_ggplot_marginal_effects.r") #plot epsilon and or survival
source("f_management_aggregate_comparison.r")
source("f_management_annual_comparison.r")
source("f_tau_sim.r")
source("f_range_sim.r")

f_cpp("v11_6") #link the library

raw_data <- "C:/NOAA/LARGE_data/chin_survival_1998_2019_with_tagged_above.csv"

sys <- Sys.time()

#Bypass models to test
bypass_mods = data.frame(t_bypass = c(0), 
                         j_bypass = c(0), 
                         l_bypass = c(0), 
                         jl_bypass = c(0), #not implemented 
                         jlt_bypass = c(0), #deprecated
                         mu_bypass = c(1)) #The full by-pass model does not converge because of the H vector

#AR models to test
AR_mods = (data.frame(t_AR = c(1), 
                     j_AR = c(2), 
                     l_AR = c(2))) 

#Random effects to test
re_mods = (data.frame(t_flag = c(1), 
                     j_flag = c(1), 
                     l_flag = c(1), 
                     jl_flag = c(0), 
                     jlt_flag = c(1)))

for(bp in 1:nrow(bypass_mods)){
  for(ar in 1:nrow(AR_mods)){
    for(re in 1:nrow(re_mods)){
      print(paste(bp,ar,re))
      #Create model combinations
      #This is super jinky. Fix your classes
      b_tmp <- as.vector(as.matrix(bypass_mods[bp,]))
      names(b_tmp) <- names(bypass_mods)
      a_tmp <- as.vector(as.matrix(AR_mods[ar,]))
      names(a_tmp) <- names(AR_mods)
      r_tmp <- as.vector(as.matrix(re_mods[re,]))
      names(r_tmp) <- names(re_mods)
      
      fit <- try(f_model(raw_file = raw_data,
                     
                     #Data trimming
                     rangeL = c(minL = 85, maxL = 125),
                     rangeJ = c(minJ = 85, maxJ = 160),
                     
                     #Number of knots
                     n_knots = 250,
                     
                     #Management actions for projecting
                     #j in days, l in mm
                     m_proj = data.frame('j' = c(-7,0,7),
                                         'l' = c(4,0,-4)),
                     
                     #Bypass offset (0 = FALSE, 1 = TRUE)
                     bypass_flags = b_tmp,
                     
                     #Random effects (0 = leave out, 1 = include)
                     re_flags = r_tmp,

                     
                     #R.E. statistical model (0 = i.i.d, 1 = AR1, 2 = R.W.) 
                     AR_flags = a_tmp,
                     
                     #TMB (random or fixed)
                     random = c("y_re", "l_re",#'z_jl', 
                                "j_re", "z_jlt"),
                     
                     H_flag = 1, #flag for anisotropy (0 = FALSE, 1 = TRUE)
                     version = "v11_6", #model version
                     compare_AIC = FALSE, #compare the AIC to previous model runs
                     getsd = FALSE, #must be turned on for marginal plots, but turned off for management plots
                     save_to_file = TRUE, #save the fit to a file
                     save_file = "bias_sim_re_y_j_l_jlt_RW_j_l_bypassmu.rData", #file name of the saved fit
                     DHARMa_sim = FALSE, #do the DHARMa simulations
                     bias_sim = TRUE, #Simulated bias of parameters
                     bias_sim_n = 50, #
                     sim_size = 1, #Sample size experiment (deprecated in cpp)
                     proj_sim = 0, #Projection simulation
                     proj_H = c(0.,0.,))) #Anisotropy simulation
    }
  }
}
  


print(Sys.time()-sys)
print(fit$opt$par)
