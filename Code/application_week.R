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

nodes <- c("d", "j")
spod <- read.csv(sprintf("../SPod/SPod_week/S08_2018-06-%d.csv",i+13), 
                 header=TRUE,  na.strings = "N/A")


spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

spodPIDs <- as.data.frame(spod[, paste(nodes, "PID..counts.", sep=".")])
names(spodPIDs) <- nodes
spodPIDs[,nodes] <- spodPIDs[,nodes]/1000
spodPIDs$time <- spod$time
spodPIDs[, nodes[1]] <- na.locf(spodPIDs[,nodes[1]])
spodPIDs[, nodes[2]] <- na.locf(spodPIDs[,nodes[2]])

spod_trends <- data.frame(time = spod$time)
window_size <- 3600
overlap <- 600
max_iter <- 5
tau <- c(0.01, 0.05, 0.1)

for (node in nodes){
  result <- get_windows_BIC(y=spodPIDs[, node], tau, k=3, window_size, overlap,
                            df_tol = 1e-9,
                            gamma = 1,
                            plot_lambda = TRUE,
                            eps_abs = 0.02, 
                            eps_rel = 1e-3,
                            solver = NULL,
                            criteria = "eBIC", 
                            max_iter = max_iter, 
                            rho = 1, 
                            update = 1)
  spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
  names(spod_trends)[(ncol(spod_trends)-length(tau)+1):ncol(spod_trends)] <-
    paste(node, tau, sep = "_")
  save(spod_trends,  spodPIDs, result,
     file = sprintf("../SPod/Results/trends_%s_2018-06-%d.RData",node,i+13))
}


