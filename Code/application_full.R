################################################################################
# Fit quantile trends on full application dataset
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
load_all("detrendr")
rm(list=ls())
spod <- read.csv("../SPod/fhrdata_2017-11-30.csv", 
                 header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

window_size <- 10000
overlap <- 1000
max_iter <- 30
tau <- c(0.05, 0.1, 0.15)
k <- 3
spod_trends <- data.frame(time = spod$time)

max_iter <- 100
result <- get_trend_windows(spodNode$pid[10001:28000], tau, 
                          lambda = 2e5,
                          k, window_size, overlap,
                          max_iter = max_iter, 
                          update = 1,
                          eps_abs = 0.005,
                          rho = 20)

result0 <- get_trend(spodNode$pid[30001:33100], tau, 
                     lambda = 1e4,
                     k)
plot(spodNode$pid[30001:33100], type="l")
lines(result0[,3])
lines(result[,3], col="red")


for (node in c("f", "g", "h")){
  pidCol <- paste(node, "SPOD.PID..V.", sep=".")
  spodNode <- spod[, c("time", pidCol)]
  names(spodNode)[2] <- c("pid")
  spodNode$pid <- spodNode$pid/1000
  result <- get_windows_BIC(spodNode$pid[30000:40000], tau, k, window_size, overlap,
                          lambdaSeq = window_size^seq(1.8, 2.3, length.out=4)[2],
                          df_tol = 1e-9,
                          gamma = 1,
                          plot_lambda = TRUE,
                          solver = NULL,
                          criteria = "eBIC", 
                          max_iter = max_iter, 
                          eps_abs = 0.01,
                          rho = 15)
  save(result, file=sprintf("../SPod/node_%s_trend.RData", node))
  spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
  names(spod_trends)[(ncol(spod_trends)-length(tau)+1):ncol(spod_trends)] <-
    paste(node, tau, sep = "_")
}
save(spod_trends, tau, file = "../SPod/spod_trends.RData")
################################################################################

