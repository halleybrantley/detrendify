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
spod <- read.csv("../SPod/fhrdata_2017-11-30.csv", 
                 header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

load("../SPod/spod_trends.RData")
nodes <- c("f", "g", "h")
spodPIDs <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")]/1000)
names(spodPIDs) <- nodes
spodPIDs$time <- spod$time

spodPIDs <- spodPIDs %>% 
  filter(time > as.POSIXct("2017-11-30 10:15:00"),
         time <= as.POSIXct("2017-11-30 10:15:00")+8000)
plot(spodPIDs$h, type="l")

tau <- c(0.05, 0.1)
k <- 3
detrendr_trends <- data.frame(time = spodPIDs$time)
qsreg_trends <- data.frame(time = spodPIDs$time)
n <- nrow(spodPIDs)
x <- seq(1, nrow(spodPIDs), 1)
for (node in c("f", "g", "h")){
  spodPIDs[, node] <- na.locf(spodPIDs[, node])
  result <- get_trend_BIC(spodPIDs[, node], tau, k, 
                            lambdaSeq = n^seq(.8, 1.5, length.out=10),
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
    fit_qsreg <- qsreg(x, spodPIDs[, node], maxit.cv = 50, 
                       alpha=tau[j], hmin = -5)
    trend[,j] <- predict(fit_qsreg)    
  }
  qsreg_trends <- cbind(qsreg_trends, as.data.frame(trend))
  names(qsreg_trends)[(ncol(qsreg_trends)-(length(tau)-1)):ncol(qsreg_trends)] <-
    paste(node, tau, sep = "_")
}

save(spodPIDs, detrendr_trends, qsreg_trends, 
     file = "../SPod/trends_short.RData")

# spodPeaks <- select(spodPIDs, -time) - select(qsreg_trends, 
#                                               contains(paste(0.1)))
# plot(spodPeaks$f, type="l")
# abline(h=0.1, col="red")
# 
# 
# thresh <- 4*apply(spodPeaks[3000:3299,], 2, sd)
# 
# lines(trend[,j], col="red")
# 
# ################################################################################
# get_spod_signal <- function(tau, spod_trends, spodPIDs, thresholds){
#   spodPeaks <- select(spodPIDs, -time) - select(spod_trends, contains(paste(tau)))
#   spodSignal <- spodPeaks
#   for (i in 1:length(thresholds)){
#     spodSignal[,i] <- as.numeric(spodPeaks[,i]>thresholds[i])
#   }
#   return(spodSignal)
# }
# 
# get_confusion <- function(spod_signal){
#   spod_signal <- na.omit(spod_signal)
#   mat_h1 <- confusionMatrix(factor(spod_signal$f[spod_signal$h==1]),
#                             factor(spod_signal$g[spod_signal$h==1]))
#   
#   mat_h0 <- confusionMatrix(factor(spod_signal$f[spod_signal$h==0]),
#                             factor(spod_signal$g[spod_signal$h==0]))
#   
#   conf_out <- cbind(mat_h0$table, mat_h1$table)
#   return(conf_out)
# }
# 
# thresh <- 0.1
# detrendr_signal <- get_spod_signal(0.1, detrendr_trends, spodPIDs, thresh)
# get_confusion(detrendr_signal)
# qsreg_signal <- get_spod_signal(0.1, qsreg_trends, spodPIDs, thresh)
# get_confusion(qsreg_signal)
# 
# latex(get_confusion(detrendr_signal),
#       file = "../Manuscript/short_confusion_detrend.tex",
#       rowlabel = "",
#       rowname = c("f = 0", "f = 1"),
#       cgroup = c("h = 0", "h = 1"),
#       colheads = rep(c("g = 0", "g = 1"),2),
#       n.cgroup = c(2,2),
#       caption = "Confusion matrices for 3 SPod nodes after baseline removal.")
# 
