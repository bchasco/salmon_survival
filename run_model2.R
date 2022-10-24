library(TMB)
library(INLA)
library(dplyr)
library(tidyr)

source("f_raw_data.r") #Collect the raw data
source("f_mesh.r") #Build the mesh
source("f_map.r") #Build the mesh
source("f_TMB_data_pars.r") #Build the data and parameters lists
source("f_model.r") #Sned everything to a TMb object and estimate
source("f_cpp.r") #Sned everything to a TMb object and estimate

f_cpp("v11_6") #link the library

fit <- f_model(raw_file = "C:/NOAA/LARGE_data/chin_survival_1998_2019_subset.csv",
                   rangeL = c(minL = 85, maxL = 125),
                   rangeJ = c(minJ = 85, maxJ = 160),
                   bypass_flags = c(mu_bypass = 0, t_bypass = 0, j_bypass = 0, 
                                    l_bypass = 0, jlt_bypass = 0),
                   H_flag = 1,
                   re_flags = c(t_flag = 1, j_flag = 1,
                                l_flag = 1, jlt_flag = 1),
                   version = "v11_6",
                   n_knots = 50,
                   random = c("y_re", "l_re", "j_re", "z_jlt"))
rm(list = ls())
