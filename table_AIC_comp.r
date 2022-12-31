load("AIC_comp.rData") #This is save output from the run_model_comb.r script
AIC_comp <- as.data.frame(AIC_comp)
names(AIC_comp) <- paste("model",1:ncol(AIC_comp))
AIC_comp <- t(t(AIC_comp)[order(AIC_comp['AIC',]),])
print("Top three models")
print(AIC_comp[,1:5])
write.csv(file="table_AIC_comp_output.csv",AIC_comp[,1:5])

