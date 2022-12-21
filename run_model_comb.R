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
bypass_mods = data.frame(t_bypass = c(1,1,0,1,1,0), 
                         j_bypass = c(0,1,0,0,1,0), 
                         l_bypass = c(0,1,0,0,1,0), 
                         jl_bypass = c(0,0,0,0,0,0), #not implemented 
                         jlt_bypass = c(0,0,0,0,0,0), #deprecated
                         mu_bypass = c(0,0,0,1,1,1)) #The full by-pass model does not converge because of the H vector

#AR models to test
AR_mods = (data.frame(t_AR = c(0,1,2), 
                     j_AR = c(0,1,2), 
                     l_AR = c(0,1,2))) 

#Random effects to test
re_mods = (data.frame(t_flag = c(0,1,1), 
                     j_flag = c(0,1,1), 
                     l_flag = c(0,1,1), 
                     jl_flag = c(0,0,0), 
                     jlt_flag = c(0,0,1)))

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
                     bypass_flags = b_tmp,#as.vector(bypass_mods[bp,]),
                     # tmp1 <- as.vector(as.matrix(bypass_mods[bp,]))
                     # names(tmp1) <- names(bypass_mods)
                     # tmp <- c(t_bypass = 1, j_bypass = 0, l_bypass = 0,
                     #                  jl_bypass = 0, #not implemented
                     #                  jlt_bypass = 0, #deprecated
                     #                  mu_bypass = 0)# #The full by-pass model does not converge because of the H vector
                     
                     
                     #Random effects (0 = leave out, 1 = include)
                     re_flags = r_tmp,#unclass(as.matrix(re_mods[re,])),
                       # c(t_flag = 1, j_flag = 1, l_flag = 1,
                       #            jl_flag = 0, jlt_flag = 1),
                     
                     
                     #R.E. statistical model (0 = i.i.d, 1 = AR1, 2 = R.W.) 
                     AR_flags = a_tmp,#unclass(as.matrix(AR_mods[ar,])),
                     # c(t_AR = 0, j_AR = 0, l_AR = 0), #AR1 with all by-pass does not converge, iid does, rw does
                     
                     #TMB (random or fixed)
                     random = c("y_re", "l_re",#'z_jl', 
                                "j_re", "z_jlt"),
                     
                     H_flag = 1, #flag for anisotropy (0 = FALSE, 1 = TRUE)
                     version = "v11_6", #model version
                     compare_AIC = TRUE, #compare the AIC to previous model runs
                     getsd = FALSE, #must be turned on for marginal plots, but turned off for management plots
                     save_to_file = FALSE, #save the fit to a file
                     save_file = ".rData", #file name of the saved fit
                     DHARMa_sim = FALSE, #do the DHARMa simulations
                     bias_sim = FALSE, #Simulated bias of parameters
                     bias_sim_n = 50, #
                     sim_size = 2, #Sample size experiment
                     proj_sim = 0, #Projection simulation
                     proj_H = c(0.,0.,))) #Anisotropy simulation
    }
  }
}
  

# source("f_create_all_plots.r")

#If you decide to change something dramatic about the model
clear_AIC_comp <- FALSE
if(clear_AIC_comp){
  AIC_comp <- list()
  save(file="AIC_comp.rData",AIC_comp)
}

source("table_AIC_comp.r")

print(Sys.time()-sys)
print(fit$opt$par)
