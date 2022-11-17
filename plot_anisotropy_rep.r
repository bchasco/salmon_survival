rep <- fit$obj$rep

Eigen = eigen(rep$H_jl)
rss = function(V) sqrt(sum(V[1]^2+V[2]^2))
Pos_Major = Eigen$vectors[,1]*Eigen$values[1] * rep$Range_raw_jl
Pos_Minor = Eigen$vectors[,2]*Eigen$values[2] * rep$Range_raw_jl
Range = 1.1 * c(-1,1) * max(abs( cbind(Pos_Major,Pos_Minor) ))

# png("plot_anisotropy.png", unit="in", height=5, width = 5, res = 600)
plot( 1, type="n", xlim=Range, ylim=c(Range[1],Range[2]*1.2), xlab="", ylab="")
shape::plotellipse( rx=rss(Pos_Major), ry=rss(Pos_Minor), angle=-1*(atan(Pos_Major[1]/Pos_Major[2])/(2*pi)*360-90), lcol="black", lty="solid")
title( "Distance at 10% correlation" )
mtext(side=1, outer=FALSE, line=2, text="Julian day")
mtext(side=2, outer=FALSE, line=2, text="Length (mm)")
# dev.off()
