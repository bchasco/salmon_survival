library(ggplot2)

#This figure is based on the model with the lowest AIC.
#It compares the aggregate estimated survival with the observed survival

# load("modSD.rData")
# load("modREP.rData")

load("sdList_v11.rData")
load("repList_v11.rData")
bm <- 3

tmp_data <- aggregate(list(s=data$surv,total=data$total),by=list(t=data$t_i+1998),sum)
tmp_data <- data.frame(val=tmp_data$s/tmp_data$total, y=tmp_data$t)
sdEst <- sdList$est[[bm]]
sdSd <- sdList$sd[[bm]]

noActions <- sdEst$proj_surv[c(1:3,17:19),9]

tmp_df <- data.frame(val = c(sdEst$proj_surv[c(1:3,17:19),1:9]),
                     sd = c(sdSd$proj_surv[c(1:3,17:19),1:9]),
                     # data = rep(tmp_data$val,times=2*9),
                     # Model = rep(c("No interaction", "Interaction"), each = length(c(mySd$sd.est$proj_surv))),
                     Action = rep(c("-7 days, +4 mm", "-7 days, -4 mm", "+7 days, -4 mm", 
                                    "+7 days, +4 mm", "+7 days, 0 mm", "-7 days, 0 mm",  
                                    "0 days, +4 mm", "0 days, -4 mm", "0 days, 0 mm"), 
                                  each=length(c(1:3,17:19))),
                     ActionOrder = as.character(rep(c(1, 3, 9 , 7, 8 , 2,  4, 6 , 5), 
                                                    each=length(c(1:3,17:19)))),
                     y = rep((1998:2016)[c(1:3,17:19)],times=9))

facet_names <- c("-7 days, +4 mm", "-7 days, -4 mm", "+7 days, -4 mm", 
                 "+7 days, +4 mm", "+7 days, 0 mm", "-7 days, 0 mm",  
                 "0 days, +4 mm", "0 days, -4 mm", "0 days, 0 mm")
names(facet_names) <- as.character(c(1, 3, 9 , 7, 8 , 2,  4, 6 , 5))

g <- ggplot(tmp_df, aes(y=(exp(val)-exp(noActions))/exp(noActions)*100, x=as.factor(y))) +
  # geom_line(size=1) +
  facet_wrap(~ActionOrder, nrow=3, labeller=labeller(ActionOrder=facet_names)) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Percent change in smolt to adult survival") +
  xlab("Migration year") + 
  geom_point() +
  # geom_point(data=tmp_data[c(1:3,17:19),],aes(x=as.factor(y), y=val),colour="red", size=4, shape=16, alpha=0.3) +
  # geom_point(aes(x=as.factor(y), y=exp(val)),colour="black", size=1, shape=15)+
  geom_errorbar(aes(ymin=(exp(val-0.67*sd)-exp(noActions))/exp(noActions)*100, 
                    ymax=(exp(val+0.67*sd)-exp(noActions))/exp(noActions)*100), width=.2,
              position=position_dodge(.9)) +
# geom_errorbar(aes(ymin=(exp(val-1.96*sd)-tmp_data[c(1:3,17:19),]$val)/tmp_data[c(1:3,17:19),]$val), 
#               ymax=(exp(val+1.96*sd)-tmp_data[c(1:3,17:19),]$val)/tmp_data[c(1:3,17:19),]$val, width=.2,
#               position=position_dodge(.9))

  geom_abline(intercept = 0, slope = 0) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

png("ggplot_management_annual_comparison.png", width=500, height = 500, res=100)

print(g)

dev.off()
