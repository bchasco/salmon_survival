par(mfrow=c(2,2))
plot(obs$x[obs$c=="length"], obs$s[obs$c=="length"]/obs$n[obs$c=="length"])
mean_jef <- 0
mean_tef <- 0
if(data$j_flag)
  mean_jef <- rep$j_re[39] *exp(opt$par["ln_sigj"])
if(data$t_flag)
  mean_tef <- rep$y_re[10] *exp(opt$par["ln_sigy"])
lines(obs$x[obs$c=="length"],plogis(opt$par["mu"]+rep$l_re*exp(opt$par["ln_sigl"])), col="black")
lines(obs$x[obs$c=="length"],plogis(opt$par["mu"]+rep$l_re*exp(opt$par["ln_sigl"]) + mean_tef + mean_jef), col="blue")
lines(obs$x[obs$c=="length"],exp(rep$s_l), col="red")

plot(obs$x[obs$c=="day"], obs$s[obs$c=="day"]/obs$n[obs$c=="day"])
mean_lef <- 0
mean_tef <- 0
if(data$l_flag)
  mean_lef <- rep$l_re[24] *exp(opt$par["ln_sigl"])
if(data$t_flag)
  mean_tef <- rep$y_re[10] *exp(opt$par["ln_sigy"])
lines(obs$x[obs$c=="day"],plogis(opt$par["mu"]+rep$j_re*exp(opt$par["ln_sigj"])), col="black")
lines(obs$x[obs$c=="day"],plogis(opt$par["mu"]+rep$j_re*exp(opt$par["ln_sigj"]) + mean_tef + mean_lef), col="blue")
lines(obs$x[obs$c=="day"],exp(rep$s_j), col="red")

plot(obs$x[obs$c=="year"], obs$s[obs$c=="year"]/obs$n[obs$c=="year"])
lines(obs$x[obs$c=="year"],plogis(opt$par["mu"]+rep$y_re*exp(opt$par["ln_sigy"])), col="black")
lines(obs$x[obs$c=="year"],plogis(opt$par["mu"]+rep$y_re*exp(opt$par["ln_sigy"]) + mean_jef + mean_lef), col="blue")
lines(obs$x[obs$c=="year"],exp(rep$s_t), col="red")

mean(rep(rep$j_re[data$j_i+1],data$total))
mean(rep(rep$l_re[data$l_i+1],data$total))
mean(rep(rep$y_re[data$t_i+1],data$total))

hist(rep(rep$j_re[data$j_i+1],data$total))
hist(rep(rep$l_re[data$l_i+1],data$total))
hist(rep(rep$y_re[data$t_i+1],data$total))
hist(rep$z_jlt/exp(opt$par["log_tau2_jl"]))


ave_l = aggregate(list(l=rep(data$l_i,data$total)+minL), by=list(t=rep(data$t_i,data$total)), mean)

load("j_ex.rdata")
load("l_ex.rdata")
load("y_ex.rdata")


# matplot(1998:2016,y_ex, type="l", lwd=2,col=1:9, lty=1:9)
#run fig2_marginalEffect.r
matplot(1998:2016,y_ex[,c("y_wozwlwj","y_wzwlwj")], type="l", col=1:4, lty=1:4)
matpoints(1998:2016,obs$s[obs$c=="year"]/obs$n[obs$c=="year"], pch=16, cex=2, col="black")
legend(12+1997,0.02,legend=names(y_ex[,c("y_wozwlwj","y_wzwlwj")]), lty=1:9, col=1:9)

matplot(85:155,j_ex[,c("j_wozwlwy","j_wzwlwy")], type="l", col=1:4, lty=1:4)
# matplot(85:155,j_ex, type="l", col=1:9, lty=1:9)
matpoints(85:155,obs$s[obs$c=="day"]/obs$n[obs$c=="day"], pch=16, cex=1.5, col="black")
legend(40+85,0.015,legend=names(j_ex[,c("j_wozwlwy","j_wzwlwy")]), lty=1:9, col=1:9)

matplot(minL:maxL,l_ex[,c("l_wozwjwy","l_wzwjwy")], type="l", col=1:4, lty=1:4)
matpoints(minL:maxL,obs$s[obs$c=="length"]/obs$n[obs$c=="length"], pch=16, cex=1.5, col="black")
legend(90,0.01,legend=names(l_ex[,c("l_wozwojwoy","l_wozwjwy","l_wzwjwy","l_wozwojwy")]), lty=1:4, col=1:4)
# legend(90,0.01,legend=c("l_re only", "l_re, j_re, y_re","all + interaction","l_re, y_re"), lty=1:4, col=1:4)

#The effect of no interaction - these assume an average year, not the years that were observed.
plot(minL:maxL,
     ylim=c(-10,13),
     (l_ex$l_wozwjwy-l_ex$l_wzwjwy)/l_ex$l_wzwjwy*100, las=1,
     ylab="Percent change in survival with no interaction term",
     xlab="Length",
     type="l", lwd=2)
abline(h=0)


s(y) = plogis(mu + re_y(y))


plot(85:155,
     (j_ex$j_wozwlwy-j_ex$j_wzwlwy)/j_ex$j_wzwlwy*100, las=1,
     ylab="Percent change in survival with no interaction term",
     xlab="Day",
     type="l", lwd=2)
abline(h=0)

plot(1998:2016,
     (y_ex$y_wozwlwj-y_ex$y_wzwlwj)/y_ex$y_wzwlwj*100, las=1,
     ylab="Percent change in survival with no interaction term",
     xlab="Year",
     type="l", lwd=2)
abline(h=0)

j_obs <- rep(data$j_i+1,data$total)
j_dist <- dnorm(1:(max(data$j_i+1)),mean(j_obs),sd(j_obs))

j_mu <- rep(data$j_i,data$total)
# j_obs_plus7 <- rep(data$j_i+1+7,data$total)
# j_mu_plus7[j_mu_plus7>max(j_obs)] <- max(j_obs)
# j_mu_minus7 <- rep(data$j_i+1-7,data$total)
# j_mu_minus7[j_mu_minus7<1] <- 1
# 
# j_timing_ex_woz <- data.frame(ex=rep(c("late","early"),each=length(j_obs)),s=c((j_ex$j_wozwlwy[j_mu_plus7]-j_ex$j_wozwlwy[j_obs])/j_ex$j_wozwlwy[j_obs],(j_ex$j_wozwlwy[j_mu_minus7]-j_ex$j_wozwlwy[j_obs])/j_ex$j_wozwlwy[j_obs]))
# j_timing_ex_woz$z <- 0
# j_timing_ex_wz <- data.frame(ex=rep(c("late","early"),each=length(j_obs)),s=c((j_ex$j_wozwlwy[j_mu_plus7]-j_ex$j_wozwlwy[j_obs])/j_ex$j_wzwlwy[j_obs],(j_ex$j_wzwlwy[j_mu_minus7]-j_ex$j_wzwlwy[j_obs])/j_ex$j_wzwlwy[j_obs]))
# j_timing_ex_wz$z <- 1

j_mu <- rep(data$j_i+1+7,data$total)
# j_mu <- rnorm(sum(data$total),mean(j_obs),sd(j_obs))
# j_obs_plus7 <- rnorm(sum(data$total),mean(j_obs)+7,sd(j_obs))
# j_mu_plus7[j_mu_plus7>max(j_obs)] <- max(j_obs)
# j_obs_minus7 <- rnorm(sum(data$total),mean(j_obs)-7,sd(j_obs))
# j_mu_minus7[j_mu_minus7<1] <- 1
# j_timing_ex_woz <- data.frame(ex=rep(c("late","early"),each=sum(data$total)),s=c((j_ex$j_wozwlwy[j_mu_plus7]-j_ex$j_wozwlwy[j_obs])/j_ex$j_wozwlwy[j_obs],(j_ex$j_wozwlwy[j_mu_minus7]-j_ex$j_wozwlwy[j_obs])/j_ex$j_wozwlwy[j_obs]))
# j_timing_ex_woz$z <- 0
# j_timing_ex_woz <- data.frame(ex=rep(c("late","early"),each=length(j_obs)),s=c((j_ex$j_wozwlwy[j_mu_plus7]-j_ex$j_wozwlwy[j_obs])/j_ex$j_wozwlwy[j_obs],(j_ex$j_wozwlwy[j_mu_minus7]-j_ex$j_wozwlwy[j_obs])/j_ex$j_wozwlwy[j_obs]))
# j_timing_ex_wz$z <- 1

j_timing_ex_woz <- data.frame(ex=rep(c("late","early"),each=length(j_mu)),s=c((j_ex$j_wozwlwy[j_mu_plus7]-j_ex$j_wozwlwy[j_mu])/j_ex$j_wozwlwy[j_mu],(j_ex$j_wozwlwy[j_mu_minus7]-j_ex$j_wozwlwy[j_mu])/j_ex$j_wozwlwy[j_mu]))
j_timing_ex_woz$z <- 0
j_timing_ex_wz <- data.frame(ex=rep(c("late","early"),each=length(j_mu)),s=c((j_ex$j_wozwlwy[j_mu_plus7]-j_ex$j_wozwlwy[j_mu])/j_ex$j_wzwlwy[j_mu],(j_ex$j_wzwlwy[j_mu_minus7]-j_ex$j_wzwlwy[j_mu])/j_ex$j_wzwlwy[j_mu]))
j_timing_ex_wz$z <- 1


j_timing_ex <- rbind(j_timing_ex_woz,
                     j_timing_ex_wz)


l_mu <- rep(data$l_i+1,data$total)
l_mu_plus10 <- rep(data$l_i+1+10,data$total)
l_mu_plus10[l_mu_plus10>max(l_mu)] <- max(l_mu)
l_mu_minus10 <- rep(data$l_i+1-10,data$total)
l_mu_minus10[l_mu_minus10<1] <- 1

l_timing_ex_woz <- data.frame(ex=rep(c("late","early"),each=length(l_mu)),s=c((l_ex$l_wozwjwy[l_mu_plus10]-l_ex$l_wozwjwy[l_mu])/l_ex$l_wozwjwy[l_mu],(l_ex$l_wozwjwy[l_mu_minus10]-l_ex$l_wozwjwy[l_mu])/l_ex$l_wozwjwy[l_mu]))
l_timing_ex_woz$z <- 0
l_timing_ex_wz <- data.frame(ex=rep(c("late","early"),each=length(l_mu)),s=c((l_ex$l_wzwjwy[l_mu_plus10]-l_ex$l_wzwjwy[l_mu])/l_ex$l_wzwjwy[l_mu],(l_ex$l_wzwjwy[l_mu_minus10]-l_ex$l_wzwjwy[l_mu])/l_ex$l_wzwjwy[l_mu]))
l_timing_ex_wz$z <- 1

l_timing_ex <- rbind(l_timing_ex_woz,
                     l_timing_ex_wz)

boxplot(j_timing_ex$s~j_timing_ex$z+j_timing_ex$ex)
boxplot(l_timing_ex$s~l_timing_ex$z+l_timing_ex$ex)


j_sim <- replicate(20,sample(rep(data$j_i,data$total),15000))
woz <- apply(j_sim,2,function(x){return(mean(j_ex$j_wozwlwy[x]))})
wz <- apply(j_sim,2,function(x){return(mean(j_ex$j_wzwlwy[x]))})
sim_diff <- data.frame(var=rep("timing",100),diff=(woz - wz)/wz)


l_sim <- replicate(20,sample(rep(data$l_i,data$total),15000))
woz <- apply(l_sim,2,function(x){return(mean(l_ex$l_wozwjwy[x]))})
wz <- apply(l_sim,2,function(x){return(mean(l_ex$l_wzwjwy[x]))})
sim_diff <- rbind(sim_diff,data.frame(var=rep("length",100),diff=(woz - wz)/wz))

boxplot(sim_diff$diff~sim_diff$var)

woz_ag_y <- aggregate(list(wt_s=nu_ex$woz*data$total, total=data$total), by=list(t=data$t_i), sum)
woz_ag_y$ex <- "woz"
wz_ag_y <- aggregate(list(wt_s=nu_ex$wz*data$total, total=data$total), by=list(t=data$t_i), sum)
wz_ag_y$ex <- "wz"
ag_y <- rbind(woz_ag_y,wz_ag_y)

library(lattice)
xyplot(data=ag_y,(wt_s/total)~t|ex)

aggregate(rep$ln_zexp_sp,by=list(data$t_i), mean)

# save(nu_i, file="nu_i.rData")
myDay <- c(120)
myLength <- c(110)
myYear <- c(1998:2016)
subDay <- rep(data$j_i+min(x$julian),2)%in%myDay
subLen <- rep(data$l_i+minL,2)%in%myLength
subYear <- rep(data$t_i+1998,2)%in%myYear
subData <- apply(cbind(subDay,subYear),1,prod)
print(sum(subData))
nu_i_ex <- data.frame(l_i=rep(data$l_i+minL,2),
                      t_i=rep(data$t_i+1998,2),
                      j_i=rep(data$j_i+min(x$julian),2),
                      pred=rep$nu_i,
                      obs = rep(data$surv/data$total,2),
                      ex=rep(c("w_oz","w_z"),each=data$n_i))

g <- ggplot(nu_i_ex[subDay==TRUE & subYear==TRUE,], aes(y=pred, x=l_i, colour=j_i, shape=ex)) +
  geom_line() +
  facet_wrap(~t_i, ncol=4) +
  geom_point(aes(x=l_i, y=obs), size=2, shape=16)+
  ylim(0,0.2)
print(g)


g <- ggplot(nu_i_ex[subLen==TRUE & subYear==TRUE,], aes(y=pred, x=j_i,  group=ex)) +
  geom_line() +
  facet_wrap(~t_i, ncol=4) +
  geom_point(aes(x=j_i, y=obs, colour="grey"), size=2, shape=16)+
  ylim(0,0.5)
print(g)

par(mfcol=c(3,2))

plot(sd_j$value[names(sd$value)=="y_re"])
lines(sd_jlt$value[names(sd$value)=="y_re"])
plot(sd_j$value[names(sd$value)=="l_re"])
lines(sd_jlt$value[names(sd$value)=="l_re"])
plot(sd_j$value[names(sd$value)=="j_re"])
lines(sd_jlt$value[names(sd$value)=="j_re"])

plot(sd_j$sd[names(sd$value)=="y_re"], main="Year")
lines(sd_jlt$sd[names(sd$value)=="y_re"])
plot(sd_j$sd[names(sd$value)=="l_re"], main="Length")
lines(sd_jlt$sd[names(sd$value)=="l_re"])
plot(sd_j$sd[names(sd$value)=="j_re"], main="Day")
lines(sd_jlt$sd[names(sd$value)=="j_re"])


par(mfcol=c(3,2))

plot(sd_j$value[names(sd$value)=="y_re"])
lines(sd_jt$value[names(sd$value)=="y_re"])
plot(sd_j$value[names(sd$value)=="l_re"])
lines(sd_jt$value[names(sd$value)=="l_re"])
plot(sd_j$value[names(sd$value)=="j_re"])
lines(sd_jt$value[names(sd$value)=="j_re"])

plot(sd_j$sd[names(sd$value)=="y_re"], main="Year")
lines(sd_jt$sd[names(sd$value)=="y_re"])
plot(sd_j$sd[names(sd$value)=="l_re"], main="Length")
lines(sd_jt$sd[names(sd$value)=="l_re"])
plot(sd_j$sd[names(sd$value)=="j_re"], main="Day")
lines(sd_jt$sd[names(sd$value)=="j_re"])
