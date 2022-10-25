load("AIC_comp.rData")
AIC_comp <- as.data.frame(AIC_comp)
names(AIC_comp) <- paste("model",1:ncol(AIC_comp))
print(t(t(AIC_comp)[order(AIC_comp['AIC',]),]))

