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
source("application_functions.R")
rm(list=ls())

spod <- read.csv("../SPod/SPod_week/SENTINEL Data_2017-03-02.csv",
                 header=TRUE,  na.strings = "N/A")
load(sprintf("../SPod/SPod_week/trends_e_2017-03-02.RData"))
load(sprintf("../SPod/SPod_week/qsreg_trends_2017-03-02.RData"))

spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 
nodes <- c("c", "e")
spodPIDs <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")])
names(spodPIDs) <- nodes
spodPIDs[,nodes] <- spodPIDs[,nodes]/1000
spodPIDs$time <- spod$time
spodPIDs[, nodes[1]] <- na.locf(spodPIDs[,nodes[1]])
spodPIDs[, nodes[2]] <- na.locf(spodPIDs[,nodes[2]])
spodPIDs <- spodPIDs[1:nrow(qsreg_trends), ]
spod_trends <- spod_trends[1:nrow(qsreg_trends), ]
tau <- c(0.01, 0.05, 0.1)

plot(spodPIDs$c, type="l")
lines(spod_trends$c_0.05, col="red")
lines(qsreg_trends$c_0.05, col="blue")

plot(spodPIDs$e, type="l")
lines(spod_trends$e_0.05, col="red")
lines(qsreg_trends$e_0.05, col="blue")

peaks_qsreg <- select(spodPIDs, -time) - 
  select(qsreg_trends, ends_with(paste(0.05)))

peaks_detrend <- select(spodPIDs, -time) - 
  select(spod_trends, ends_with(paste(0.05)))

plot(c~e, peaks_qsreg)
plot(c~e, peaks_detrend)
cor(peaks_detrend[,c("c", "e")], method="spearman")
cor(peaks_qsreg[,c("c", "e")], method="spearman")


plot(peaks_detrend$e, type="l")
lines(peaks_detrend$c, col="red")

plot(peaks_qsreg$e, type="l")
lines(peaks_qsreg$c, col="red")

node <- "c"



getCutoff(peaks_detrend, "e")


signal_detrend <- peaks_detrend
for (node in nodes){
  signal_detrend[,node] <- as.numeric(signal_detrend[,node]>
                                        getCutoff(signal_detrend, node))
}

signal_qsreg <- peaks_qsreg
for (node in nodes){
  signal_qsreg[,node] <- as.numeric(signal_qsreg[,node]>
                                      getCutoff(signal_qsreg, node))
}

get_VI(signal_detrend, nodes)
get_VI(signal_qsreg, nodes)

VI_d <- VI_q <- pos <- c()
for (j in 1:4){
  ind_start <- (j-1)*21600 + 1
  ind_end <- min(nrow(spodPIDs), j*21600)
  pos[j] <- sum(apply(cbind(signal_detrend[ind_start:ind_end,], 
                            signal_qsreg[ind_start:ind_end,]), 1, 
                      function(x) any(x>0)))
  VI_d[j] <- get_VI(signal_detrend[ind_start:ind_end,], nodes)
  VI_q[j] <- get_VI(signal_qsreg[ind_start:ind_end,], nodes)
}
cbind(VI_d, VI_q)

