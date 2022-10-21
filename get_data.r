data <- read.csv("PIT_survival.csv")
data <- read.csv("D:/CRB Survival analysis/raw_spring_chin_data.csv")
data$y_i <- data$year
data$d_i <- data$doy
data$l_i <- data$length

ag_data <- aggregate(list(count=data$survived), #This is a 0 or 1
                     by=list(y=data$y_i, d=data$d_i, l=data$l_i, b=data$barge_binary, temp=data$daily_mean_Temp, tdg=data$daily_mean_tdg, spill=data$AvgTotalSpill, surv=data$survived),
                     FUN=function(x){return(length(x))})

library(reshape2)
d_ag <- dcast(ag_data, y + d + b + l + temp + tdg + spill ~ surv, value.var="count")
d_ag$`0`[is.na(d_ag$`0`)] <- 0
d_ag$`1`[is.na(d_ag$`1`)] <- 0
d_ag$total <- d_ag$`0` + d_ag$`1`
names(d_ag)[8:9] <- c('dead','surv') 
save(d_ag, file="d_ag_env.Rdata")
d_ag_env<-d_ag


####Create environemental indices 
flow <- read.csv("D:/CRB Survival analysis/environmental/LGR_flow.csv")
temp <- read.csv("D:/CRB Survival analysis/environmental/LGR_temp_tdg.csv")

flow$Date2 <- as.Date(flow$Date,format='%m/%d/%Y')
flow$doy <- lubridate::yday(flow$Date2)
flow$year <- lubridate::yday(flow$Date2)

library(plyr)

temp_daily<-ddply(temp, .(day), summarize, dm_Temp = mean(temp_c), dm_tdg=mean(tdg))
temp_daily$Date2 <- as.Date(temp_daily$day,format='%m/%d/%Y')
temp_daily$doy <- lubridate::yday(temp_daily$Date2)

head(temp_daily)

enviro <- merge(x= temp_daily, y= flow, by.x= 'day',by.y='Date', all.x= T)
enviro<-subset(enviro, select=-c(Date2.y, doy.y,Date2.x, doy.x))
enviro$Date2 <- as.Date(enviro$day,format='%m/%d/%Y')
enviro$doy <- lubridate::yday(enviro$Date2)
enviro$year <- lubridate::year(enviro$Date2)

library(zoo)
enviro$discharge_20mean<-rollmean(enviro$AvgTotalDischarge, 20, na.pad = TRUE, align = "right")
enviro$temp_20mean<-rollmean(enviro$dm_Temp, 20, na.pad = TRUE, align = "right")
enviro$tdg_20mean<-rollmean(enviro$dm_tdg, 20, na.pad = TRUE, align = "right")


