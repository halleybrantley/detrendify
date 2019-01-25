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

window_size <- 5000
overlap <- 500
max_iter <- 20
tau <- c(0.05, 0.1, 0.15)
k <- 3
spod_trends <- data.frame(time = spod$time)


for (node in c("f", "g", "h")){
  pidCol <- paste(node, "SPOD.PID..V.", sep=".")
  spodNode <- spod[, c("time", pidCol)]
  names(spodNode)[2] <- c("pid")
  spodNode$pid <- spodNode$pid/1000
  result <- get_windows_BIC(spodNode$pid, tau, k, window_size, overlap,
                          lambdaSeq = window_size^seq(0.8, 1.5, length.out=10),
                          df_tol = 1e-9,
                          gamma = 1,
                          plot_lambda = FALSE,
                          solver = NULL,
                          criteria = "eBIC", 
                          max_iter = max_iter, 
                          eps_abs = 0.01, 
                          eps_rel = 0.002, 
                          update = 1)
  save(result, file=sprintf("../SPod/node_%s_trend.RData", node))
  spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
  names(spod_trends)[(ncol(spod_trends)-length(tau)+1):ncol(spod_trends)] <-
    paste(node, tau, sep = "_")
}
save(spod_trends, tau, file = "../SPod/spod_trends.RData")
################################################################################

