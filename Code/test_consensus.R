library(devtools)
load_all("detrendr")
rm(list=ls())
n <- 1000
x <- seq(1, n, 1)
y <- sin(x*2*pi/n) + rnorm(n, 0, .4)
lambda <- n
k <- 3
tau <- c(0.05, 0.5)
overlap <- 20
rho <- 1
window_size <- 100
result <- consensus_ADMM(y, tau, lambda, k, rho, window_size, overlap, 300, 0.05)

y_n <- length(y)
window_size <- floor((y_n+overlap)/2)
plot(result$phiBar[,1]~x, type="l")
points(result$phi1[,1]~x[1:window_size])
points(result$phi2[,1]~x[(y_n-window_size+1):y_n], col="blue")
plot(y~x, type="l", col="grey")
lines(result$theta[,1], col="red")
lines(result$theta[,2], col="blue")
plot(result$primal_norm)
plot(result$dual_norm)
result$iter

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

plot(pid~time, spodNode[35000:36500, ], type="l")

y <- spodNode[34801:36600, "pid"]
x <- spodNode[34801:36600, "time"]

k <- 3
tau <- c(0.1, 0.5)
overlap <- 300
rho <- 1
window_size <- 900
lambda <- 10*window_size

result <- consensus_ADMM(y, tau, lambda, k, rho, window_size, overlap, 200, 0.1)
y_n <- length(y)

plot(y, type="l")
lines(result$theta[,1], col="red")
lines(result$theta[,2], col="blue")
abline(v=window_size)
abline(v=window_size-overlap, col="purple")
abline(v=2*window_size-overlap, col="purple")
abline(v=2*window_size-2*overlap, col="darkgreen")
abline(v=3*window_size-2*overlap, col="darkgreen")

plot(result$primal_norm[20:40])
plot(result$dual_norm[20:80])
result$iter

plot(y~x, type="l")
lines(theta[,1]~x, col="red")
lines(theta[,2]~x, col="blue")

plot(phiBar[,2]~x, type="l")
points(phi1[,2]~x[1:window_size])
points(phi2[,2]~x[(y_n-window_size+1):y_n], col="blue")
