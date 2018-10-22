# Spingarn Experiments

library(Rglpk)
library(plyr)
library(devtools)
library(microbenchmark)
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
spodNode$pid <- scale(spodNode$pid, center = TRUE)

timeBase <- "10 sec"
timeBreaks <- seq(round(min(spodNode$time), "min"), 
                  round(max(spodNode$time), "min"), timeBase)
spodNode$timeCut <- cut(spodNode$time, timeBreaks)
spod10 <- ddply(subset(spodNode, !is.na(timeCut)), .(timeCut), 
                summarize, 
                pid = mean(pid, na.rm=TRUE))


k <- 3
tau <- .2
lambda <- 600
y <- spod10$pid[1:200]
D <- as.matrix(get_Dk(length(y), k))

microbenchmark(
theta_gurobi <- gurobi_lp(y, tau, lambda, D),
theta_lpgl <- lpglpk_trendfilter(y, tau, lambda, D))

plot(y, type="l")
lines(theta_gurobi, col="red")
lines(theta_lpgl, col="blue")

theta_lpgl <- list()

n <- 1000
i <- 0
for (lambda in c(1, 10, 100, 1000)){
  i <- i + 1
  y <- spod10$pid[1:n]
  D <- as.matrix(get_Dk(length(y), k))
  theta_lpgl[[i]] <- lpglpk_trendfilter(y, tau, lambda, D)
}


plot(y, type="l", col="grey", main = "LP solver")
lines(theta_lpgl[[1]], col = "red")
lines(theta_lpgl[[2]], col = "orange")
lines(theta_lpgl[[3]], col = "blue")
lines(theta_lpgl[[4]], col = "green")

# df <- data.frame(n = c(200, 400, 800), 
#                  cl = NA, 
#                  normD = NA)
# 
# i <- 0
# for (n in df$n){
#   i <- i + 1
#   y <- spod10$pid[1:n]
#   D <- as.matrix(get_Dk(length(y), k))
#   df$cl[i] <- check_loss(y-theta_lpgl[[i]], tau)
#   df$normD[i] <- norm(D%*%theta_lpgl[[i]], "1")
# 
# }



step_size <- 1
theta0 <- y
eta0 <- D%*%theta0 
M <- Cholesky(Diagonal(n) + Matrix::crossprod(D))

x0 <- c(theta0, as.numeric(eta0))
y0 <- double(length(x0))

admm <- ADMM(x0, y0, prox_f, prox_g, step_size = .0001,
             max_iter = 5000, eps_abs=5e-6, eps_rel = 1e-6,
             1000, y, D, M, tau)

x0 <- admm$x
y0 <- admm$y
plot(y, type="l")
lines(admm$x[1:n], col="blue", type="l")
lines(theta_lpgl[[4]], col="red")

# Plot residuals 
plot(admm$dual_norm[200:1000], type="l")
plot(admm$primal_norm[200:1000], type="l")
lines(admm$eps_dual, col="red")
lines(admm$eps_pri, col="red")


# spign_const <- spingarn_multi_step_R(theta0, eta0, y, D, M, lambda, 
#                                      tau, 1, numberIter=10000, 
#                                      function(x, y){x})
# plot(y, type="l")
# lines(admm$x[1:n], col="blue", type="l")
# lines(theta_lpgl[[4]], col="red")
# lines(spign_const$theta)

#FASTA
fasta_fit <- fasta(f, gradf, g, proxg, x0, tau1 = 1, max_iters = 10000, 
                   w = 10, backtrack = TRUE, recordIterates = TRUE, 
                   stepsizeShrink = 0.5, eps_n = 1e-15, 
                   nu=.0001, y=y, D=D, tau = tau, lambda = n)

plot(y, type="l")
lines(theta_lpgl[[i]], col="red", type="l")
lines(fasta_fit$x[1:n], col="blue")

################################################################################
# Linearized ADMM
x0 <- theta0 
y0 <- double(length(D%*%x0))

pr_f1 <- function(x, step, y, tau, lambda){
  prox_f1(x, y, tau, step)
}

pr_f2 <- function(x, step, y, tau, lambda){
  prox_f2(x, lambda, step)
} 

ladmm <- LADMM(x0, y0, D, pr_f1, pr_f2, step_f = 1, 
               step_g = norm(D, "f")^2,
               max_iter = 20000,
               y, tau, lambda)

x0 <- ladmm$x
y0 <- ladmm$y

plot(y, type="l", col="grey")
lines(theta_lpgl[[3]], col="red")
lines(ladmm$x, col="cyan")


