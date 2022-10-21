library(ggplot2)

load("sdList_v11.rData")
load("repList_v11.rData")
bm <- 3

tmp_data <- aggregate(list(s=data$surv,total=data$total),by=list(t=data$t_i+1998),sum)
tmp_data <- data.frame(val=tmp_data$s/tmp_data$total, y=tmp_data$t)
sdEst <- sdList$est[[bm]]
sdSd <- sdList$sd[[bm]]



tmp_df <- data.frame(val = c(sdEst$proj_surv_ag),
                     sd = c(sdSd$proj_surv_ag),
                     Action = c("-7 days, +4 mm", 
                                "-7 days, -4 mm", 
                                "+7 days, -4 mm", 
                                "+7 days, +4 mm", 
                                "+7 days, 0 mm", 
                                "-7 days, 0 mm",  
                                "0 days, +4 mm", 
                                "0 days, -4 mm", 
                                "0 days, 0 mm"),
                     ActionOrder = c(1, 3, 9 , 7, 8 , 2,  4, 6 , 5),
                     y = NA)
Action = c("-7 days, +4 mm", "-7 days, -4 mm",
           "+7 days, -4 mm", "+7 days, +4 mm",
           "+7 days, 0 mm", "-7 days, 0 mm",
           "0 days, +4 mm", "0 days, -4 mm",
           "0 days, 0 mm")
ActionOrder = c(1, 3, 9 , 7, 8 , 2,  4, 6 , 5)


g <- ggplot(tmp_df, aes(y=exp(val), x=as.factor(ActionOrder))) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Smolt to adult survival") +
  xlab("Management action") + 
  geom_errorbar(aes(ymin=exp(val-1.96*sd), ymax=exp(val+1.96*sd)), width=.2,
                position=position_dodge(.9)) +
  geom_point(aes(x=as.factor(ActionOrder), y=exp(val)),colour="grey", size=4, shape=15)+
  geom_hline(yintercept = sum(data$surv)/sum(data$total), colour="red", size=3, alpha=0.3) +
  theme(axis.text.x = element_text(angle = 90)) +
  # geom_point(data=tmp_data[c(1:3,17:19),],aes(x=as.factor(y), y=val),colour="red", size=4, shape=16, alpha=0.3)
  scale_x_discrete(labels= Action[order(ActionOrder)])


png("ggplot_management_aggregate_comparison.png", width=500, height = 500, res=100)

print(g)

dev.off()
