################################################################################
# Fit quantile trends on full application dataset
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
library(zoo)
load_all("detrendr")
rm(list=ls())
spod <- read.csv("../SPod/raw/SENTINEL Data_2017-04-13.csv", 
                 header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

surge0 <- which(spod$d.SPOD.Sonic.Voltage > 14)
surge <- c()
for (i in seq(-70, 70, 1)){
  surge <- c(surge, surge0+i)
}
surge <- unique(surge)

nodes <- c("c", "d", "e")
spodPIDs <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")])
names(spodPIDs) <- nodes
spodPIDs$time <- spod$time
spodPIDs$d[surge] <- NA
for (node in nodes){
  spodPIDs[,node] <- spodPIDs[,node]/1000 
}
save(spodPIDs, file = "../SPod/spodPIDs.RData")

window_size <- 3600
overlap <- 600
max_iter <- 10
tau <- c(0.01, 0.05, 0.1)
k <- 3
spod_trends <- data.frame(time = spod$time)

for (node in c("c", "d", "e")){
  missID <- which(is.na(spodPIDs[, node]))
  spodPIDs[,node] <- na.approx(spodPIDs[,node], na.rm=FALSE)
  spodPIDs[missID, node] <- spodPIDs[missID, node]  + rnorm(length(missID), 0, .002)
  spodNode <- spodPIDs[, c("time", node)]
  names(spodNode)[2] <- c("pid")
  result <- get_windows_BIC(spodNode$pid, tau, k, window_size, overlap,
                          lambdaSeq = exp(c(4, 5, 6, 7, 8, 9, 10, 11, 12, 13)),
                          df_tol = 1e-9,
                          gamma = 1,
                          plot_lambda = TRUE,
                          solver = NULL,
                          criteria = "eBIC", 
                          max_iter = max_iter, 
                          rho = 1, 
                          update = 2)
  save(result, file=sprintf("../SPod/node_%s_trend.RData", node))
  spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
  names(spod_trends)[(ncol(spod_trends)-length(tau)+1):ncol(spod_trends)] <-
    paste(node, tau, sep = "_")
}
save(spod_trends, tau, file = "../SPod/Results/trends_e_2017-04-13.RData")
################################################################################

