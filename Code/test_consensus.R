library(devtools)
library(doParallel)
build("detrendr")
install("detrendr")
load_all("detrendr")
library(detrendr)
library(tidyverse)
library(profvis)

rm(list=ls())
# Application
dataDir <- "~/Desktop/EPA/SPod_Data/TestRange_Dec2017"
datafiles <- dir(dataDir, pattern=".csv", full.names=TRUE)
spod <- read.csv(datafiles[3], header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

node <- "h"
pidCol <- paste(node, "SPOD.PID..V.", sep=".")
spodNode <- spod[, c("time", pidCol)]
names(spodNode)[2] <- c("pid")
spodNode <- subset(spodNode, !is.na(pid))
spodNode$pid <- as.numeric(scale(spodNode$pid, center = TRUE))

plot(pid~time, spodNode[35000:37000, ], type="l")

y <- spodNode[35001:37000, "pid"]
x <- (seq(1, length(y), 1)-1)/2000 #spodNode[33801:37000, "time"]
n <- length(y)
k <- 3
tau <- c(0.05, 0.5)
overlap <- 60
rho <- 1
window_size <- 1200
lambda <- 5*length(y)
max_iter <- 100

result <- consensus_ADMM(y, tau, lambda, k, rho, window_size, overlap, 
                                     max_iter, eps = 2e-4, update = 1)

result0 <- gurobi_trend(y, tau, lambda, k)

finit <- aresult(x,(2000-100+1),0.04,hry1(100,2000,x),hry2(100,tau,2000,y))  #Yu's method
fit_Oh <- qreq1d.sreg(x,y,finit,p=tau,cutoff=0.001) 

y_n <- length(y)

plot(y, type="l", col="grey")
lines(result0[,1], col="red")
lines(fit_Oh, col="blue")

lines(result$theta[,1], col="blue")
abline(v=window_size, col="darkgreen")
abline(v=window_size-overlap, col="purple")

n_windows <- ceiling(length(y)/(window_size-overlap))

for (i in 2:n_windows){
  abline(v=i*window_size-(i-1)*overlap, col="purple")
  abline(v=i*window_size-i*overlap, col="darkgreen")
}

dev.copy(png, "consensus_fig.png", width = 800, height = 400)
dev.off()

theta.df <- get_windows(y, x, k, tau, lambda, window_size, overlap)

