################################################################################
# Fit smoothing spline quantile trends (qsreg) on week of data
# Halley Brantley
################################################################################
library(tidyverse)
library(fields)
library(zoo)
rm(list=ls())
tau <- c(0.01, 0.05, 0.1)
load("../SPod/spodPIDs.RData")

nodes <- c("c", "d", "e")
for (node in nodes){
  missInd <- which(is.na(spodPIDs[,node]))
  spodPIDs[,node] <- na.locf(spodPIDs[,node]) 
  spodPIDs[missInd, node] <- spodPIDs[missInd, node] + 
    rnorm(length(missInd), 0, .001)
}  
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
save(qsreg_trends, tau, nodes, spodPIDs,
    file = "../SPod/Results/qsreg_trends_2017-04-13.RData")
