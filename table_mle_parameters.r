load('fit_model_55.rData')

df <- data.frame(mle = fit$opt$SD$value[grep("tab",names(fit$opt$SD$value))],
           sd = fit$opt$SD$sd[grep("tab",names(fit$opt$SD$value))],
           mu = 0,
           lwr = 0,
           upr = 0)

df$mu[1:9] <- exp(df$mle[1:9])
df$mu[10:15] <- plogis(df$mle[10:15])
df$lwr[1:9] <- exp(df$mle[1:9] - 1.28  * df$sd[1:9])
df$lwr[10:15] <- plogis(df$mle[10:15] - 1.28  * df$sd[10:15])
df$upr[1:9] <- exp(df$mle[1:9] + 1.28  * df$sd[1:9])
df$upr[10:15] <- plogis(df$mle[10:15] + 1.28  * df$sd[10:15])

df <- df[df$sd>0,]

print(df)
