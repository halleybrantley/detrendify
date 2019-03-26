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

i = 0 
spod <- read.csv(sprintf("../SPod/SPod_week/SENTINEL Data_2017-03-0%d.csv",i), 
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

spod_trends <- data.frame(time = spod$time)
window_size <- 3600
overlap <- 600
max_iter <- 5
tau <- c(0.1, 0.15)


for (node in nodes){
  result <- get_windows_BIC(y=spodPIDs[,node], tau, k=3, window_size, overlap,
                            lambdaSeq = exp(seq(5,13,1)),
                            df_tol = 1e-9,
                            gamma = 1,
                            plot_lambda = TRUE,
                            solver = NULL,
                            criteria = "eBIC", 
                            max_iter = max_iter, 
                            rho = 1, 
                            update = 1)
  spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
  names(spod_trends)[(ncol(spod_trends)-length(tau)+1):ncol(spod_trends)] <-
    paste(node, tau, sep = "_")
  save(spod_trends, spodPIDs, result, 
       file = sprintf("../SPod/SPod_week/trends_2017-03-0%d.RData",i))
}


