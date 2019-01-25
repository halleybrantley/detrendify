################################################################################
# Fit quantile trends on shorter window
# Halley Brantley
################################################################################
library(fields)
library(tidyverse)
library(devtools)
library(Hmisc)
library(caret)
library(zoo)
load_all("detrendr")
rm(list=ls())
spod <- read.csv("../SPod/fhrdata_2017-11-30.csv", 
                 header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

load("../SPod/spod_trends.RData")
nodes <- c("f", "g", "h")
spodPID <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")]/1000)
names(spodPID) <- nodes
spodPID$time <- spod$time

spodPIDs <- spodPID %>% 
  filter(time > as.POSIXct("2017-11-30 9:30:00"),
         time <= as.POSIXct("2017-11-30 9:30:00")+5000)
plot(spodPIDs$h, type="l")

tau <- c(0.05, 0.1, 0.15)
k <- 3
detrendr_trends <- data.frame(time = spodPIDs$time)
qsreg_trends <- data.frame(time = spodPIDs$time)
n <- nrow(spodPIDs)
x <- seq(1, nrow(spodPIDs), 1)
for (node in c("f", "g", "h")){
  spodPIDs[, node] <- na.locf(spodPIDs[, node])
  result <- get_trend_BIC(spodPIDs[, node], tau, k, 
                            lambdaSeq = n^seq(1.1, 1.8, length.out=10),
                            df_tol = 1e-9,
                            gamma = 1,
                            plot_lambda = TRUE,
                            solver = "gurobi",
                            criteria = "eBIC")
  detrendr_trends <- cbind(detrendr_trends, as.data.frame(result$theta))
  names(detrendr_trends)[(ncol(detrendr_trends)-(length(tau)-1)):
                           ncol(detrendr_trends)] <-
    paste(node, tau, sep = "_")
  
  trend <- matrix(NA, n, length(tau))
  for (j in 1:length(tau)){
    fit_qsreg <- qsreg(x, spodPIDs[, node], maxit.cv = 50, 
                       alpha=tau[j], hmin = -6)
    trend[,j] <- predict(fit_qsreg)    
  }
  qsreg_trends <- cbind(qsreg_trends, as.data.frame(trend))
  names(qsreg_trends)[(ncol(qsreg_trends)-(length(tau)-1)):ncol(qsreg_trends)] <-
    paste(node, tau, sep = "_")
}

save(spodPIDs, detrendr_trends, qsreg_trends, tau,
     file = "../SPod/trends_short.RData")

