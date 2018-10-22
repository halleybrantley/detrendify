# Spingarn Experiments

library(Rglpk)
load_all("detrendr")

rm(list=ls())

dataDir <- "~/Desktop/EPA/SPod_Data/TestRange_Dec2017"
datafiles <- dir(dataDir, pattern=".csv", full.names=TRUE)
spod <- read.csv(datafiles[3], header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

node <- "f"
pidCol <- paste(node, "SPOD.PID..V.", sep=".")
spodNode <- spod[, c("time", pidCol)]
names(spodNode)[2] <- c("pid")
spodNode <- subset(spodNode, !is.na(pid))
spodNode <- spodNode[35001:45000,]
plot(pid~time, spodNode, type="l", main = "raw")

y <- spodNode$pid/1000
n <- length(y)
k <- 1
tau <- .05
D <- as.matrix(get_Dk(length(y), 1))
lambda <- 1
theta_lpgl <- lpglpk_trendfilter(y, tau, lambda, D)
step_size = 1e-6

theta0 <- rep(quantile(y, .5), length(y))
eta0 <- D%*%theta0 + rnorm(nrow(D), sd = .01)
M <- Cholesky(Diagonal(n) + Matrix::crossprod(D))

spign_const <- spingarn_multi_step_R(theta0, eta0, y, D, M, lambda, 
                               tau, 1, numberIter=10000, 
                               function(x, y){x})
  
spign_C2 <- spingarn_multi_step_R(theta0, eta0, y, D, M, lambda, 
                                  tau, step_size, numberIter=10000, 
                                  function(x, y){x + 50/(y+50)})

spign_C3 <- spingarn_multi_step_R(theta0, eta0, y, D, M, lambda, 
                                  tau, step_size, numberIter=10000, 
                                  function(x, y){x + exp(-y/200)})

x0 <- c(theta0, as.numeric(eta0))
y0 <- double(length(x0))

admm2 <- ADMM2(x0, y0, f, g, lambda, 1, y, D, M, tau, max_iter = 100)
x0 <- admm2$x
y0 <- admm2$y
plot(x0[1:n], type="l")
lines(theta_lpgl)

k <- 3
D <- as.matrix(get_Dk(length(y), k))
lambda <- 2
theta_lpgl <- lpglpk_trendfilter(y, tau, lambda, D)
M <- Cholesky(Diagonal(n) + Matrix::crossprod(D))
admm2 <- ADMM2(x0, y0, f, g, lambda, 1, y, D, M, tau, max_iter = 100)
x0 <- admm2$x
y0 <- admm2$y
plot(x0[1:n], type="l")
lines(theta_lpgl)

plot(log(spign_const$diff_theta), type="l", ylim=c(-19, -5))
lines(log(spign_C2$diff_theta), col="orange")
lines(log(spign_C3$diff_theta), col="blue")

plot(log(spign_const$steps), type="l", ylim=c(-5, 0))
lines(log(spign_C2$steps), col="orange")
lines(log(spign_C3$steps), col="blue")

plot(y, type="l", col="grey")
lines(theta_lpgl, lty = 2, type="l")
lines(spign_const$theta)
lines(spign_C2$theta, col="orange")
lines(spign_C3$theta, col="blue")
lines(admm_x0[1:n], col="purple")
lines(admm2$x[1:n], col="blue")
lines(y, col="grey")

mean(abs(spign_const$theta - theta_lpgl))
mean(abs(spign_C2$theta - theta_lpgl))
mean(abs(spign_C3$theta - theta_lpgl))
mean(abs(admm2$[1:n] - theta_lpgl))

theta0 <- spign_C2$theta
eta0 <- spign_C2$eta


plot(D%*%theta_lpgl)
plot(D%*%theta0)
