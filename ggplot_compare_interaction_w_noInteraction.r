load("simout.Rdata")

load("repList.rData")
load("sdList.rData")
bm <- 3
myList <- repList[[bm]]

Action = c("-7 days, +4 mm", "-7 days, -4 mm",
           "+7 days, -4 mm", "+7 days, +4 mm",
           "+7 days, 0 mm", "-7 days, 0 mm",
           "0 days, +4 mm", "0 days, -4 mm",
           "0 days, 0 mm")
ActionOrder = c(1, 3, 9 , 7, 8 , 2,  4, 6 , 5)

facet_names <- c("-7 days, +4 mm", "-7 days, -4 mm", "+7 days, -4 mm", 
                 "+7 days, +4 mm", "+7 days, 0 mm", "-7 days, 0 mm",  
                 "0 days, +4 mm", "0 days, -4 mm", "0 days, 0 mm")
names(facet_names) <- as.character(c(1, 3, 9 , 7, 8 , 2,  4, 6 , 5))

df_x <- data.frame(nox = c(repList[[3]]$proj_surv_nox),
                   x = c(repList[[3]]$proj_surv),
                   Action = rep(Action,19),
                   ActionOrder = rep(ActionOrder,19))
                   
p <- ggplot(df_x, aes(x=plogis(nox), y = plogis(x))) + 
  geom_point()+
  facet_wrap(~ActionOrder, nrow=3, labeller=labeller(ActionOrder=facet_names))+
  geom_abline(
    slope = 1,
    intercept = 0,
  )+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  ylab("With interaction") +
  xlab("Without interaction") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

  

# png("ggplot_compare_interaction_w_noInteraction.png", height=500, width = 500, res=100)

print(p)

# dev.off()

  
