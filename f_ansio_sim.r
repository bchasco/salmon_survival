o.H <- fit$obj$env$last.par[names(fit$obj$env$last.par) == 'ln_H_input_jl'] 
proj_H <- matrix(c(o.H,c(0,0),c(0,0)),3,2, byrow = TRUE)


set.seed(100)
proj_sim_n <- 20
proj_sim <- array(NA, c(proj_sim_n,3,c(dim(fit$obj$rep$proj_y))))
for(j in 1:nrow(proj_H)){
  for(i in 1:proj_sim_n){
    fit$obj$env$data$proj_sim <- 1
    fit$obj$env$last.par[which(names(fit$obj$env$last.par) == 'ln_H_input_jl')] <- proj_H[j,]
    proj_sim[i,j,,,] <- fit$obj$simulate(complete=TRUE)$proj_y
  }
}

x0.5 <- round(sum(na.omit(c(proj_sim[,1,1:2,3,])<0))/length(na.omit(c(proj_sim[,1,1:2,3,])))*100,1)
df <- data.frame(val = c(proj_sim[,1,1:2,3,]),
                 var=rep(paste0("anisotropic (", x0.5,"%)"),length(c(proj_sim[,1,1:2,3,]))))

x1 <- round(sum(na.omit(c(proj_sim[,2,1:2,3,])<0))/length(na.omit(c(proj_sim[,2,1:2,3,])))*100,1)
df <- rbind(df,
            data.frame(val = c(proj_sim[,2,1:2,3,]),
                       var=rep(paste0("isotropic (", x1,"%)"),length(c(proj_sim[,2,1:2,3,])))))


png(file="f_H_sim.png")
g <- ggplot(df, aes(x=val, color=as.factor(var))) +
  geom_histogram(aes(fill=as.factor(var)), alpha=0.5, position="identity") +
  xlab("Percent change in annual survival") +
  ylab("") +
  xlim(-50,100) +
  geom_vline(xintercept = 0) +
  guides(fill=guide_legend(title="Interaction\ncorrelation"))

print(g)
dev.off()

fit$obj$env$last.par[names(fit$obj$env$last.par) == 'ln_H_input_jl'] <- o.H
