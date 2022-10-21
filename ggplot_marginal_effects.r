library(ggplot2)
library(ggthemes) # Load
library(cowplot)

#Year
myObjects <- c("sd_none", "sd_jlt")
myRE <- c("mar_y","mar_l", "mar_j")
xx <- list(y=1998:2016,
           l = minL:maxL,
           j = 1:length(parameters$j_re)+min(x$julian)-1)
xNames <- c(y="Year",
           l = "Length (mm)",
           j = "Migration past Lower Granite Dam")
xintercept <- c(mean(rep(data$t_i,data$total))+1998,
               mean(rep(data$l_i+minL,data$total)),
               mean(rep(data$j_i+min(x$julian),data$total)))


icnt <- 1
for(i in myRE){
  df1 <- data.frame(x = rep(xx[[icnt]],2),
                    val = c(get(myObjects[1])$value[names(sd$value)==i],get(myObjects[2])$value[names(sd$value)==i]),
                    sd = c(get(myObjects[1])$sd[names(sd$value)==i], get(myObjects[2])$sd[names(sd$value)==i]),
                    Model = c(rep("No interaction",length(get(myObjects[1])$sd[names(sd$value)==i])), rep("Interaction",length(get(myObjects[2])$sd[names(sd$value)==i]))))

  
  p <- ggplot(df1, aes(x=x, y = plogis(val), colour=Model))+
    geom_line(size=1.1) +
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Survival") +
      theme(legend.position = "none")+
    xlab(xNames[icnt]) +
    # geom_hline(yintercept = 
    #              mean(rep(data$surv/data$total,data$total)))+
    # geom_vline(aes(xintercept =2007)+
    geom_ribbon(aes(ymin=plogis(val-1.96*sd), ymax=plogis(val+1.96*sd), fill=Model), linetype = 0, alpha=0.25)

  if(icnt ==3){
    p <- ggplot(df1, aes(x = as.Date(x, origin = as.Date("2018-01-01")), y = plogis(val), colour=Model))+
      geom_line(size=1.1) +
      theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
      ylab("Survival") +
      theme(legend.position = "none")+
      xlab(xNames[icnt]) +
      # geom_hline(yintercept = 
      #              mean(rep(data$surv/data$total,data$total)))+
      # geom_vline(xintercept = 
      #              as.Date(xintercept[icnt], origin = as.Date("2018-01-01")))+
      geom_ribbon(aes(ymin=plogis(val-1.96*sd), ymax=plogis(val+1.96*sd), fill=Model), linetype = 0, alpha=0.25)

  }
  assign(i,p)
  icnt <- icnt + 1
}

df2 <- data.frame(var = rep(c(rep("Length", length(parameters$l_re)),rep("Day", length(parameters$j_re)),rep("Year",length(parameters$y_re))),2),
                  Model = c(rep("No interaction",length(sd$value[grep("mar",names(sd$value))])), rep("Interaction",length(sd$value[grep("mar",names(sd$value))]))),
                  val = c(sd_none$sd[grep("mar",names(sd$value))], sd_jlt$sd[grep("mar",names(sd$value))]))
p <- ggplot(df2, aes(y=val, x=var, fill=Model), alpha = 0.25) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("Standard error") +
  xlab("Random effect") +
  geom_violin(trim=FALSE, alpha=0.4) 




plot_grid(get(myRE[1]), get(myRE[2]), get(myRE[3]), p, 
          labels = c("A\n106 mm\nMay 3rd", "B\n2007\nMay 3rd", "C\n2007\n106 mm", "D"),
          ncol = 2, nrow = 2)
