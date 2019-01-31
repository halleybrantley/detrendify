################################################################################
# Fit quantile trends on full application dataset
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
load_all("detrendr")
rm(list=ls())
spod <- read.csv("../SPod/SENTINEL Data_2017-04-13.csv", 
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
  spodPIDs[,node] <- (spodPIDs[,node] - min(spodPIDs[,node], na.rm=T))*5/
    (max(spodPIDs[,node], na.rm=T)-min(spodPIDs[,node], na.rm=T))
}
save(spodPIDs, file = "../SPod/spodPIDs.RData")

window_size <- 5000
overlap <- 1000
max_iter <- 30
tau <- c(0.15)
k <- 3
spod_trends <- data.frame(time = spod$time)


result <- get_trend_windows(spodNode$pid[10001:19000], tau, 
                          lambda = exp(17),
                          k, window_size, overlap,
                          max_iter = max_iter, 
                          update = 1,
                          rho = 1)

result0 <- get_trend(spodNode$pid[30001:33100], tau, 
                     lambda = 1e4,
                     k)
plot(spodNode$pid[20001:29000], type="l")
lines(result$trend[,1], col="red")
plot(result$trend[,1], col="red", type="l")


for (node in c("c", "d", "e")){
  spodNode <- spodPIDs[, c("time", node)]
  names(spodNode)[2] <- c("pid")
  result <- get_windows_BIC(spodNode$pid[20001:29000], tau, k, window_size, overlap,
                          lambdaSeq = exp(seq(12,19,1)),
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
save(spod_trends, tau, file = "../SPod/spod_trends.RData")
################################################################################

