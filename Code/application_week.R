################################################################################
# Fit quantile trends on week of data
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
library(fields)
library(gurobi)
load_all("detrendr")
rm(list=ls())

i = 2
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

qsreg_trends <- data.frame(time = spodPIDs$time, c_0.1 = NA, 
                        c_0.15 = NA,  e_0.1 = NA, e_0.15 = NA)
spod_trends <- data.frame(time = spod$time)
window_size <- 10000
overlap <- 1000
max_iter <- 5
tau <- c(0.1, 0.15)
k <- 3


tau <- c(0.1, 0.15)
t1 <- Sys.time()
for (j in 1:12){
  ind_start <- (j-1)*7200 + 1
  ind_end <- min(nrow(spodPIDs), j*7200)
  x <- seq(ind_start, ind_end, 1)
  n <- length(x)
  trends <- data.frame(time = spodPIDs$time[ind_start:ind_end])
  for (node in nodes){
    trend <- matrix(NA, n, length(tau))
    for (k in 1:length(tau)){
      fit_qsreg <- qsreg(x, spodPIDs[ind_start:ind_end,node], 
                         maxit.cv = 50, 
                       alpha=tau[k], hmin = -6, hmax = NA)
      trend[,k] <- predict(fit_qsreg)   
    }
    trends <- cbind(trends, as.data.frame(trend))
    names(trends)[(ncol(trends)-(length(tau)-1)):ncol(trends)] <-
      paste(node, tau, sep = "_")
  }
  qsreg_trends[ind_start:ind_end, ] <- trends
  
  for (node in nodes){
    result <- get_windows_BIC(spodPIDs[,node], tau, k, window_size, overlap,
                              lambdaSeq = exp(seq(12,19,1)),
                              df_tol = 1e-9,
                              gamma = 1,
                              plot_lambda = TRUE,
                              solver = NULL,
                              criteria = "eBIC", 
                              max_iter = max_iter, 
                              rho = 1, 
                              update = 2)
    spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
    names(spod_trends)[(ncol(spod_trends)-length(tau)+1):ncol(spod_trends)] <-
      paste(node, tau, sep = "_")
  }
}

save(spod_trends, qsreg_trends, spodPIDs, file = )
