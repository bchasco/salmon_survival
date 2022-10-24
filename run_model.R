library(TMB)
library(INLA)
library(dplyr)
library(tidyr)

source("f_raw_data.r")
source("f_mesh.r")
source("f_TMB_data_pars.r")

raw_file =  "C:/NOAA/LARGE_data/chin_survival_1998_2019_subset.csv"
rangeL <- c(minL = 85, maxL = 125)
rangeJ <- c(minJ = 85, maxJ = 160)

bypass_flags <- c(mu_bypass = 0, t_bypass = 0, j_bypass = 0, 
                 l_bypass = 0, jlt_bypass = 0)
re_flags <- c(t_flag = 1, j_flag = 1,
             l_flag = 1, jlt_flag = 1)
H_flag <- 1


df <- f_raw_data(file = raw_file,
                  rangeL = rangeL,
                  rangeJ = rangeJ)

mesh <- f_mesh(data = df,
               n_knots = n_knots)

TMB_list <- f_TMB_data_pars(df = df,
                            mesh = mesh,
                            by_pass = 1,
                            b_flags = bypass_flags,
                            re_flags = re_flags)
version <- "v11_6"
source("build_cpp_v1.r")


random <- c("y_re", "l_re", "j_re", "z_jlt")

map <- f_map(TMB_list)

version <- "v11_6"

obj <- MakeADFun(TMB_list$data
                 ,silent = FALSE
                 , TMB_list$parameters
                 , random=random
                 , map = map
                 , DLL=paste0("spde_aniso_wt_",version)
                 , sdreport = FALSE)

obj$mesh <- mesh

opt <- TMBhelper::fit_tmb( obj,
                           loopnum = 1,
                           newtonsteps = 1,
                           quiet = TRUE,
                           getsd = FALSE) # where Obj is a compiled TMB object

obj$rep <- obj$report()

object.size(obj)
