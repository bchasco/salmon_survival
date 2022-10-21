# library(RandomFields)
library(INLA)
library(raster)
library(RANN)

df <- data.frame(x$length,
                 x$julian,
                 x$year,
                 x$survived,
                 x$above)

names(df) <- c('len'
               ,'j'
               ,'y'
               ,'s'
               ,'a'
)

# df <- df[df$j>60 & df$j<=160,] #Get rid of some of the dates

#Aggregate data to improve speed
ag_s <- aggregate(list(s=df$s), 
                  by=list(len=df$len,j=df$j,y=df$y,a=df$a),
                  sum)
ag_t <- aggregate(list(t=df$s), 
                  by=list(len=df$len,j=df$j,y=df$y,a=df$a),
                  function(x){return(length(x))})
df2 <- data.frame(ag_s,ag_t$t)
names(df2)[ncol(df2)] <- c('ns')


#Raw locations
loc_jl <- cbind(df2$j,df2$l)
#map locations to knot centers
if( n_knots < nrow(loc_jl) ){
  knots_xy_jl = kmeans( x=loc_jl, centers=n_knots)
  # Modify data
  loc_xy_jl = knots_xy_jl$centers
  s_i_jl = knots_xy_jl$cluster
}

# Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
mesh_jl = inla.mesh.create( loc_xy_jl, refine=TRUE, extend=1 )
inla_spde_jl = inla.spde2.matern(mesh = mesh_jl)



# ---------- Start code that prepare objects for anisotropy. Should not be modified when you make your own example
Dset = 1:2
# Triangle info
TV = mesh_jl$graph$tv           # Triangle to vertex indexing
V0 = mesh_jl$loc[TV[,1],Dset]   # V = vertices for each triangle
V1 = mesh_jl$loc[TV[,2],Dset]
V2 = mesh_jl$loc[TV[,3],Dset]
E0 = V2 - V1                      # E = edge for each triangle
E1 = V0 - V2
E2 = V1 - V0  
# Calculate Areas
TmpFn = function(Vec1, Vec2) abs(det( rbind(Vec1, Vec2) ))
Tri_Area = rep(NA, nrow(E0))
for(i in 1:length(Tri_Area)) Tri_Area[i] = TmpFn( E0[i,],E1[i,] )/2   # T = area of each triangle
# ---------- End code that prepare objects for anisotropy. 

spde_jl <- list(
  "n_s"      = inla_spde_jl$n.spde,
  "n_tri"    = nrow(TV),
  "Tri_Area" = Tri_Area,
  "E0"       = E0,
  "E1"       = E1,
  "E2"       = E2,
  "TV"       = TV - 1,
  "G0"       = inla_spde_jl$param.inla$M0,
  "G0_inv"   = as(diag(1/diag(inla_spde_jl$param.inla$M0)), "dgTMatrix"))

n_x_jl <- mesh_jl$n
n_i <- nrow(df2)

