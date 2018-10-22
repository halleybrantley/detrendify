# ADMM Experiments

library(Rglpk)
load_all("../Code/detrendr")

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
spodNode <- spodNode[35001:35500,]
plot(pid~time, spodNode, type="l", main = "raw")

y <- spodNode$pid/1000
n <- length(y)
k <- 2
tau <- .05
D <- as.matrix(get_Dk(length(y), k))
lambda <- .001
theta_lpgl <- lpglpk_trendfilter(y, tau, lambda, D)

step_size = .1

theta0 <- rep(quantile(y, .5), length(y))
eta0 <- D%*%theta0 + rnorm(nrow(D), sd = .1)
x0 <- c(theta0, as.numeric(eta0))
M <- Cholesky(Diagonal(n) + Matrix::crossprod(D))

x0 <- ADMM(x0, f, g, lambda, step_size, y, D, M, tau, max_iter = 15000)

plot(y, type="l")
lines(theta_lpgl, col="red")
lines(x0[1:n], col="blue")
