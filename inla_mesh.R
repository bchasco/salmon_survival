library(INLA)
library(TMB)


#ag_chin <- read.csv("surv.csv", header=TRUE)
#Create the data 
rm(list=ls())

try(dyn.unload("inla_mesh"))
compile("inla_mesh.cpp")
dyn.load("inla_mesh")

load("d_ag_env.Rdata")
j_exp <- 0
l_exp <- 0

min_d <- 70 #lower re range
min_data_d <- 70 #data limit
max_d <- 200
max_data_d <- 200
min_l <- 65
min_data_l <- 65
max_l <- 135
max_data_l <- 125
dBin <- 10
lBin <- 10
t_ag <- d_ag[d_ag$l>=min_data_l & d_ag$l<=max_data_l,]
t_ag <- t_ag[t_ag$d>=min_data_d & t_ag$d<=max_data_d,] #
t_ag$d <- ceiling((t_ag$d - min_data_d)/dBin)
t_ag$l <- ceiling((t_ag$l - min_data_l)/lBin)

n_j <- length(min_data_d:max_data_d)
n_l <- length(min_data_l:max_l)

st <- 1
s <- 1 
re_y <- 1
# re_j <- 0
# re_l <- 0

b <- c(1)

t_ag_surv <- aggregate(list(dead=t_ag$dead, surv=t_ag$surv, total=t_ag$total),
                       by = list(y=t_ag$y, 
                                 d=t_ag$d, 
                                 b=t_ag$b, 
                                 l=t_ag$l), 
                       sum)

t_ag_env <- aggregate(list(temp=t_ag$temp-mean(t_ag$temp),tdg=t_ag$tdg-mean(t_ag$tdg)), 
                      by = list(y=t_ag$y, 
                                d=t_ag$d, 
                                b=t_ag$b, 
                                l=t_ag$l), 
                      mean)

t_ag <- cbind(t_ag_surv,t_ag_env[,c("temp","tdg")])
t_ag <- t_ag[t_ag$b%in%b,] 
t_ag_tmp <- aggregate(list(total=t_ag$total,surv=t_ag$surv),by=list(b=t_ag$b,y=t_ag$y),sum)

total_t <- matrix(t_ag_tmp$total,length(b),length(unique(t_ag$y)))
total_s <- matrix(t_ag_tmp$surv,length(b),length(unique(t_ag$y)))

it <- Sys.time()

# lattice <- inla.mesh.lattice(
#   seq(1,ceiling((max_d-min_data_d)/dBin),1),
#   seq(1,ceiling((max_l-min_data_l)/lBin),1),
# )
# mesh = inla.mesh.create(lattice=lattice,
#   extend = FALSE)
loc <- cbind(t_ag$d,t_ag$l)
mesh = inla.mesh.create(loc=loc,
                        refine=TRUE,
                        extend=-1)
inla_spde = inla.spde2.matern(mesh = mesh)

# plot(mesh)
# points(t_ag$d,t_ag$l, pch=16, cex=1.8, col="grey")
# points(loc, pch=16, cex=1.8, col="grey")


x_s <- NA
icnt <- 1
for(i in 1:nrow(t_ag)){
  x_s[icnt] <- (1:dim(mesh$loc)[1])[mesh$loc[,1]==t_ag$d[i] & mesh$loc[,2]==t_ag$l[i]]
  icnt <- icnt + 1
}

n_x <- mesh$n
n_i <- nrow(t_ag)
# n_l <- length(unique(t_ag$l - min(t_ag$l)))
# n_j <- length(unique(t_ag$d - min(t_ag$d)))
data <- list(n_i = n_i
             ,n_b = length(b)
             ,b_i = t_ag$b
             ,n_x = n_x
             ,n_t = length(unique(t_ag$y))
             ,x_s = x_s - 1#mesh$idx$loc-1
             ,proj_x_s = x_s - 1#mesh$idx$loc-1
             ,t_i = t_ag$y - min(t_ag$y)#Use all of the dataa
             ,l_i = t_ag$l - min(t_ag$l)
             ,j_i = t_ag$d - min(t_ag$d)
             ,total = t_ag$total
             ,total_t = total_t
             ,total_s = total_s
             ,surv = t_ag$surv
             ,s = s
             ,st = st
             ,re_y = re_y
             ,G0 = inla_spde$param.inla$M0
             ,G1 = inla_spde$param.inla$M1
             ,G2 = inla_spde$param.inla$M2
            
)
if(length(b)==1){
  data$b_i <- data$b_i*0
}

parameters <- list(
  mu = -4
  ,log_tau_O = rep(0,data$n_b)      # log-inverse SD of Omega
  ,log_tau_E = 0      # log-inverse SD of Omega
  ,log_kappa = 0      # Controls range of spatial variation
  ,frho = 0.
  ,frho_q = 0.
  ,Omega_input = matrix(0,data$n_x,data$n_b)   # Spatial variation in carrying capacity
  ,Epsilon_input = matrix(0,data$n_x,data$n_t)   # Spatial variation in carrying capacity
  ,eps_y = matrix(0,data$n_b,data$n_t)
  ,frho_y = 0.
  ,frho_y2 = 0.
  ,ln_sigma_y = rep(0,data$n_b)
)

random <- c("Omega_input", "Epsilon_input", "eps_y")

source("create_inla_map.r")

j_names <- unique(ceiling((min_data_d:max_d-min_d)/dBin))*dBin + min_d
l_names <- unique(ceiling((min_data_l:max_l-min_l)/lBin))*lBin + min_l
y_names <- min(t_ag$y):max(t_ag$y)

file <- paste("inla_",mod)

if(!length(grep(file,list.files()))){
  obj <- MakeADFun(data=data
                   ,parameters = parameters
                   ,random=random
                   ,map=myMap
                   ,DLL="inla_mesh") #Create the TMB object
  
  opt <- nlminb(obj$par, obj$fn, obj$gr) #Using the non-linear optimizer for the TMB object
  rep <- obj$report()
  sd <- sdreport(obj)
  sd.est <- as.list(sd,"Estimate", report=TRUE)
  sd.sd <- as.list(sd,"Std. Error", report=TRUE)
  # opt <- TMBhelper::fit_tmb(obj = obj
  #                           ,getsd = TRUE
  #                           ,newtonsteps = 1 #We hsould make this at least 2
  # )
  # save(file=paste0(file,".rData"),
  #      opt,
  #      obj,
  #      rep,
  #      mesh,
  #      lattice,
  #      y_names,
  #      l_names,
  #      j_names)
} else {
  load(file=paste0(file,".rData"))
}
  

print(Sys.time()-it)

source("ggplot_proj_s.r")
source("ggplot_dayXyear.r")
source("ggplot_st.r")
