################################################################################
# Fit quantile trends on week of data
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
library(fields)
library(gurobi)
library(zoo)
load_all("detrendr")
rm(list=ls())

for (d in 2:8){
spod <- read.csv(sprintf("../SPod/SPod_week/SENTINEL Data_2017-03-0%d.csv",d), 
                 header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 
nodes <- c("c", "e")
spodPIDs <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")])
names(spodPIDs) <- nodes
spodPIDs$c <- spodPIDs$c/1000
spodPIDs$e <- spodPIDs$e/1000
spodPIDs$time <- spod$time
spodPIDs$c <- na.locf(spodPIDs$c)
spodPIDs$e <- na.locf(spodPIDs$e)

qsreg_trends <- data.frame(time = spodPIDs$time, c_0.1 = NA, 
                           c_0.15 = NA,  e_0.1 = NA, e_0.15 = NA)
tau <- c(0.1, 0.15)

for (j in 1:12){
  ind_start <- (j-1)*7200 + 1
  ind_end <- min(nrow(spodPIDs), j*7200)
  x <- seq(ind_start, ind_end, 1)
  n <- length(x)
  trends <- data.frame(time = spodPIDs$time[ind_start:ind_end])
  for (node in nodes){
    trend <- matrix(NA, n, length(tau))
    for (i in 1:length(tau)){
      fit_qsreg <- qsreg(x, spodPIDs[ind_start:ind_end,node], 
                         maxit.cv = 50, 
                         alpha=tau[i], hmin = -12, hmax = NA)
      trend[,i] <- predict(fit_qsreg)   
    }
    trends <- cbind(trends, as.data.frame(trend))
    names(trends)[(ncol(trends)-(length(tau)-1)):ncol(trends)] <-
      paste(node, tau, sep = "_")
  }
  qsreg_trends[ind_start:ind_end, ] <- trends
 }
 save(qsreg_trends,
     file = sprintf("../SPod/SPod_week/qsreg_trends_2017-03-0%d.RData",d))
}