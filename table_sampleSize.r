myTab <- aggregate(list(Survivors=data$surv, Total = data$total), by=list(Year=data$t_i+1998),sum)
write.csv(myTab, file="table_sampleSize.csv", quote = FALSE, row.names = FALSE)
