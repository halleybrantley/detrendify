library(devtools)
library(Hmisc)
library(fields)
library(gurobi)
library(zoo)
load_all("detrendr")
rm(list=ls())
spodAll <- {}
for (i in 2:8){
  spod <- read.csv(sprintf("../SPod/SPod_week/SENTINEL Data_2017-03-0%d.csv",i), 
                 header=TRUE,  na.strings = "N/A")
  spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                   format= "%m/%d/%Y %H:%M:%S")) 
  nodes <- c("c", "e")
  spodPIDs <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")])
  names(spodPIDs) <- nodes
  spodAll <- bind_rows(spodAll, spodPIDs)
}

par(mfrow=c(2,1))
plot(spodAll$c, type="l")
plot(spodAll$e, type="l")

min(spodAll$c, na.rm=TRUE)
min(spodAll$e, na.rm=TRUE)
