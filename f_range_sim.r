# Type Range_raw_jl = sqrt(8.0) / exp( log_kappa_jl );
o.kappa <- fit$obj$env$last.par['log_kappa_jl'] 
ranges <- sqrt(8)/exp( o.kappa )
ranges <- ranges*c(0.5,1,2)
proj_kappa <- log(sqrt(8)/ranges)

set.seed(100)
proj_sim_n <- 20
proj_sim <- array(NA, c(proj_sim_n,3,c(dim(fit$obj$rep$proj_y))))
jcnt <- 1
for(j in proj_kappa){
  for(i in 1:proj_sim_n){
    fit$obj$env$data$proj_sim <- 1
    fit$obj$env$last.par[which(names(fit$obj$env$last.par) == 'log_kappa_jl')] <- j
    proj_sim[i,jcnt,,,] <- fit$obj$simulate(complete=TRUE)$proj_y
  }
  jcnt <- jcnt + 1
}

x0.5 <- round(sum(na.omit(c(proj_sim[,1,1:2,3,])<0))/length(na.omit(c(proj_sim[,1,1:2,3,])))*100,1)
df <- data.frame(val = c(proj_sim[,1,1:2,3,]),
                 var=rep(paste0("0.5x (",x0.5,"%)"),length(c(proj_sim[,1,1:2,3,]))))

x1 <- round(sum(na.omit(c(proj_sim[,2,1:2,3,])<0))/length(na.omit(c(proj_sim[,2,1:2,3,])))*100,1)
df <- rbind(df,
            data.frame(val = c(proj_sim[,2,1:2,3,]),
                       var=rep(paste0("1x (",x1,"%)"),length(c(proj_sim[,2,1:2,3,])))))

x2 <- round(sum(na.omit(c(proj_sim[,3,1:2,3,])<0))/length(na.omit(c(proj_sim[,3,1:2,3,])))*100,1)
df <- rbind(df,
            data.frame(val = c(proj_sim[,3,1:2,3,]),
                       var=rep(paste0("2x (",x2,"%)"),length(c(proj_sim[,3,1:2,3,])))))

png(file="f_kappa_sim.png")
g <- ggplot(df, aes(x=val, color=as.factor(var))) +
  geom_histogram(aes(fill=as.factor(var)), alpha=0.5, position="identity") +
  xlab("Percent change in annual survival") + 
  ylab("") +
  xlim(-50,100) +
  geom_vline(xintercept = 0) +
  guides(fill=guide_legend(title="Interaction\nrange"))

print(g)
dev.off()

fit$obj$env$last.par['log_kappa_jl'] <- o.kappa
