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

spod <- read.csv("../SPod/SPod_week/SENTINEL Data_2017-04-13.csv",
                 header=TRUE,  na.strings = "N/A")

spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 
nodes <- c("c", "e")
spodPIDs <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")])
names(spodPIDs) <- nodes
spodPIDs[,nodes] <- spodPIDs[,nodes]/1000
spodPIDs$time <- spod$time
spodPIDs[, nodes[1]] <- na.locf(spodPIDs[,nodes[1]])
spodPIDs[, nodes[2]] <- na.locf(spodPIDs[,nodes[2]])

plot(c~e, spodPIDs[1:7200,])
tau <- c(0.01, 0.05, 0.1)
j <- 3
ind_start <- (j-1)*7200 + 1
ind_end <- min(nrow(spodPIDs), j*7200)
x <- seq(ind_start, ind_end, 1)
n <- length(x)

qsreg_trends <- data.frame(time = spodPIDs$time)

for (j in 1:12){
  ind_start <- (j-1)*7200 + 1
  ind_end <- min(nrow(spodPIDs), j*7200)
  x <- seq(ind_start, ind_end, 1)
  n <- length(x)
  trends <- data.frame(time = spodPIDs$time[ind_start:ind_end])
  for (node in nodes){
    trend <- matrix(NA, n, length(tau))
    for (i in 1:length(tau)){
      fit_qsreg <- qsreg(x, spodPIDs[ind_start:ind_end,node], 
                         maxit.cv = 50, 
                         alpha=tau[i], hmin = -12, hmax = NA)
      trend[,i] <- predict(fit_qsreg)   
    }
    trends <- cbind(trends, as.data.frame(trend))
    names(trends)[(ncol(trends)-(length(tau)-1)):ncol(trends)] <-
      paste(node, tau, sep = "_")
  }
  if (j == 1){
    qsreg_trends[,names(trends)[2:ncol(trends)]] <- NA
  }
  qsreg_trends[ind_start:ind_end, ] <- trends
  print(j)
}


qsreg_trends <- qsreg_trends[1:ind_end,]
spod_trends <- data.frame(time = spod$time[1:ind_end])
window_size <- 3600
overlap <- 600
max_iter <- 5

for (node in nodes){
  result <- get_windows_BIC(y=spodPIDs[1:ind_end, node], tau, k=3, 
                            window_size, overlap,
                            lambdaSeq = c(30000, 35000, 40000, 45000),
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
}

save(spod_trends, qsreg_trends,spodPIDs, file="../SPod/trends_4-13.Rdata")
plot(c~e, spodPIDs[1:ind_end,])

plot(spodPIDs$c[1:ind_end], type="l")
lines(spod_trends$c_0.05, col="red")
lines(qsreg_trends$c_0.05, col="blue")

plot(spodPIDs$e[1:ind_end], type="l")
lines(spod_trends$e_0.05, col="red")
lines(qsreg_trends$e_0.05, col="blue")

peaks_qsreg <- select(spodPIDs[1:ind_end,], -time) - 
  select(qsreg_trends[1:ind_end,], ends_with(paste(0.05)))

peaks_detrend <- select(spodPIDs[1:ind_end,], -time) - 
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


getThresh <- function(x){mean(x, na.rm=T)+4*sd(x, na.rm=T)}

getCutoff <- function(peaks, node){
  thresh0 <- getThresh(peaks[,node])
  print(thresh0)
  peaks[which(peaks[,node] > thresh0), node] <- NA
  thresh <- getThresh(peaks[,node])
  while(thresh0 - thresh > 0.001){
    thresh0 <- thresh
    print(thresh0)
    peaks[which(peaks[,node] > thresh), node] <- NA
    thresh <- getThresh(peaks[,node])
  }
  thresh
}

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

VI_d <- VI_q <- c()
for (j in 1:10){
  ind_start <- (j-1)*3600 + 1
  ind_end <- min(nrow(spodPIDs), j*3600)
  signal_detrend <- get_spod_signal(0.05, spod_trends[ind_start:ind_end,], 
                                    spodPIDs[ind_start:ind_end,], crit = 5)
  signal_qsreg <- get_spod_signal(0.05, qsreg_trends[ind_start:ind_end,], 
                                  spodPIDs[ind_start:ind_end,], crit = 5)
  pos[j] <- sum(apply(cbind(signal_detrend, signal_qsreg), 1, 
                      function(x) any(x>0)))
  VI_d[j] <- get_VI(signal_detrend, nodes)
  VI_q[j] <- get_VI(signal_qsreg, nodes)
}
cbind(VI_d, VI_q)
VI_q
