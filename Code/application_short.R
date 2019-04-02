################################################################################
# Fit quantile trends on 2 hour window from 4-13-17
# Halley Brantley
################################################################################
library(fields)
library(tidyverse)
library(devtools)
library(Hmisc)
library(zoo)
load_all("detrendr")
rm(list=ls())

load("../SPod/spodPIDs.RData")
nodes <- c("c", "d", "e")
spodPID <- spodPIDs
spodPIDs <- spodPID %>% 
  filter(time > as.POSIXct("2017-04-13 13:10:00"),
         time <= as.POSIXct("2017-04-13 13:10:00")+7200)

plot(spodPIDs$d, type="l")
tau <- c(0.01, 0.05, 0.1)
k <- 3

detrendr_trends <- data.frame(time = spodPIDs$time)
qsreg_trends <- data.frame(time = spodPIDs$time)
n <- nrow(spodPIDs)
x <- seq(1, nrow(spodPIDs), 1)
for (node in nodes){
  missID <- which(is.na(spodPIDs[, node]))
  spodPIDs[,node] <- na.approx(spodPIDs[,node], na.rm=FALSE)
  spodPIDs[missID, node] <- spodPIDs[missID, node]  +
    rnorm(length(missID), 0, .002)
  result <- get_trend_BIC(spodPIDs[, node], tau, k,
                            df_tol = 1e-9,
                            lambdaSeq= exp(seq(5, 10, 1)),
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
    fit_qsreg <- qsreg(x, spodPIDs[,node], maxit.cv = 50,
                       alpha=tau[j], hmin = -16, hmax = NA)
    trend[,j] <- predict(fit_qsreg)
  }
  qsreg_trends <- cbind(qsreg_trends, as.data.frame(trend))
  names(qsreg_trends)[(ncol(qsreg_trends)-(length(tau)-1)):ncol(qsreg_trends)] <-
    paste(node, tau, sep = "_")
  spodPIDs[missID, node] <- NA
}

save(spodPIDs, detrendr_trends, qsreg_trends, tau,
     file = "../SPod/trends_short.RData")

load("../SPod/trends_short.RData")
