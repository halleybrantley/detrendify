# FASTA Experiments

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
k <- 3
tau <- .05
D <- as.matrix(get_Dk(length(y), k))
lambda <- 1
theta_lpgl <- lpglpk_trendfilter(y, tau, lambda, D)
nu <- .03 #n/lambda/(2)
step_size = 1


theta0 <- rep(quantile(y, .5), length(y))
eta0 <- D%*%theta0 + rnorm(nrow(D), sd = .01)
x0 <- c(theta0, as.numeric(eta0))

fasta_fit <- fasta(f, gradf, g, proxg, x0, tau1 = step_size, max_iters = 80000, 
                   w = 10, backtrack = TRUE, recordIterates = TRUE, 
                   stepsizeShrink = 0.5, eps_n = 1e-15, 
                   nu=nu, y=y, D=D, tau = tau, lambda = lambda)
plot(fasta_fit$objective[1000:length(fasta_fit$objective)], type="l")
theta0 <- fasta_fit$x[1:n]
eta0 <- fasta_fit$x[(n+1):length(x0)]
x0 <- fasta_fit$x

plot(y, type="l")
lines(theta_lpgl, col="red")
lines(theta0, col="blue")
