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
# spod <- read.csv("../SPod/fhrdata_2017-11-30.csv", 
#                  header=TRUE,  na.strings = "N/A")
# spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
#                                  format= "%m/%d/%Y %H:%M:%S")) 
# spodPID <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")]/1000)
# names(spodPID) <- nodes
# spodPID$time <- spod$time

load("../SPod/spodPIDs.RData")
nodes <- c("c", "d", "e")
spodPID <- spodPIDs
spodPID[, nodes] <- as.data.frame(scale(spodPID[, nodes]))
spodPIDs <- spodPID %>% 
  filter(time > as.POSIXct("2017-04-13 13:20:00"),
         time <= as.POSIXct("2017-04-13 13:20:00")+6000)

plot(spodPIDs$c, type="l")
tau <- c(0.5, 0.1, 0.15)
k <- 3

detrendr_trends <- data.frame(time = spodPIDs$time)
qsreg_trends <- data.frame(time = spodPIDs$time)
n <- nrow(spodPIDs)
x <- seq(1, nrow(spodPIDs), 1)
for (node in nodes){
  missID <- which(is.na(spodPIDs[, node]))
  spodPIDs[,node] <- na.locf(spodPIDs[,node])
  spodPIDs[missID, node] <- spodPIDs[missID, node]  + rnorm(length(missID), 0, .001)
  result <- get_trend_BIC(spodPIDs[, node], tau, k, 
                            lambdaSeq = exp(seq(11,18,1)),
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
    fit_qsreg <- qsreg(x, spodPIDs[,node], maxit.cv = 50, 
                       alpha=tau[j], hmin = -6, hmax = NA)
    trend[,j] <- predict(fit_qsreg)    
  }
  qsreg_trends <- cbind(qsreg_trends, as.data.frame(trend))
  names(qsreg_trends)[(ncol(qsreg_trends)-(length(tau)-1)):ncol(qsreg_trends)] <-
    paste(node, tau, sep = "_")
}

save(spodPIDs, detrendr_trends, qsreg_trends, tau,
     file = "../SPod/trends_short.RData")

