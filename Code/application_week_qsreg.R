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
tau <- c(0.01, 0.05, 0.1)

for (d in 2){
  if (d == 1){
    nodes <- c("c", "d", "e")
    spod <- read.csv("../SPod/SPod_week/SENTINEL Data_2017-04-13.csv", 
                     header=TRUE,  na.strings = "N/A")
  } else{
    nodes <- c("c", "g")
    spod <- read.csv(sprintf("../SPod/SPod_week/S08_2018-06-%d.csv",d+13), 
                     header=TRUE,  na.strings = "N/A")
  }

  spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                   format= "%m/%d/%Y %H:%M:%S")) 
  
  spodPIDs <- as.data.frame(spod[, paste(nodes, "PID..counts.", sep=".")])
  names(spodPIDs) <- nodes
  spodPIDs[,nodes] <- spodPIDs[,nodes]/1000
  spodPIDs$time <- spod$time
  spodPIDs[, nodes[1]] <- na.locf(spodPIDs[,nodes[1]])
  spodPIDs[, nodes[2]] <- na.locf(spodPIDs[,nodes[2]])
  
  qsreg_trends <- data.frame(time = spodPIDs$time)

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
    if (j == 1){
      qsreg_trends[,names(trends)[2:ncol(trends)]] <- NA
    }
    qsreg_trends[ind_start:ind_end, ] <- trends
    print(j)
  }
  if (d == 1){
    save(qsreg_trends, tau, nodes, spodPIDs,
        file = "../SPod/SPod_week/qsreg_trends_2017-04-13.RData")
  } else {
    save(qsreg_trends, tau, nodes, spodPIDs,
      file = sprintf("../SPod/SPod_week/qsreg_trends_2018-06-%d.RData",d+13))
  }
}
