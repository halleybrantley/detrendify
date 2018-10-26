# Spingarn Experiments

library(Rglpk)
library(gurobi)
library(plyr)
library(devtools)
library(microbenchmark)
library(trustOptim)
load_all("detrendr")

rm(list=ls())

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

timeBase <- "5 sec"
timeBreaks <- seq(round(min(spodNode$time), "min"), 
                  round(max(spodNode$time), "min"), timeBase)
spodNode$timeCut <- cut(spodNode$time, timeBreaks)
spod_agg <- ddply(subset(spodNode, !is.na(timeCut)), .(timeCut), 
                summarize, 
                pid = mean(pid, na.rm=TRUE))
spod_agg$time <- as.POSIXct(as.factor(spod_agg$timeCut))

k <- 3
tau <- c(0.05, 0.1, 0.5)
df <- spod_agg
y <- df$pid[1:4500]
n <- length(y)
lambdaSeq <- c(n/10, n/5, n/2, n, 2*n, 5*n, 10*n)

loss_cv <- lambda_cv(y, tau, k, lambdaSeq,
          numFolds = 10, parallel = TRUE)
plot(rowSums(loss_cv$loss))
loss_cv$lambda

theta_df <- get_windows(df$pid, df$time, k, tau, loss_cv$lambda, 
                        length(y), 0)

plot.df <- left_join(theta_df, spod_agg) 

ggplot(plot.df, aes(x = time, y = pid)) +
    geom_line(alpha = 0.2) +
    geom_line(aes(y=theta, col = factor(window), linetype=tau)) +
    theme_bw()


df$cor <- df$pid - theta_df[which(theta_df$tau== "tau_0.05"), "theta"]

# node_f <- df
# node_g <- df
ggplot(df[7000:8000,], aes(x=time, y=cor)) + geom_line() + 
  theme_bw() + 
  geom_line(data=node_f[7000:8000,], col="blue", alpha = .5)+ 
  geom_line(data=node_g[7000:8000,], col="red", alpha = .5)


#####################################################################
tau <- 0.05
lambda <- 10*length(y)
phi0 <- y- quantile(y, tau)
phiHat <- lbfgs(obj, grad_obj, vars = phi0,  y=y, tau=tau, D=D, 
                lambda = lambda, invisible = 1)
D <- get_Dk(length(y), k)
thetaHat <- y-phiHat$par
plot(y, type="l", col="grey")
lines(thetaHat)

phiHat <- trust.optim(phi0, obj, grad_obj, hess_obj, 
                      method = "Sparse", 
                      control = list(maxit = 50, 
                                     report.level = 4, 
                                     start.trust.radius = 0.04, 
                                     trust.iter = 5000), 
                      y=y, tau=tau, 
                      D=D, lambda=lambda)

thetaHat <- y-phiHat$solution
phi0 <- phiHat$solution
plot(y, type="l", col="grey")
lines(thetaHat)

theta2 <- gurobi_trend(y, tau, lambda, k)
lines(theta2, col="red")

phi <- phi0[1:10]
y <- y[1:10]
D <- get_Dk(length(y), k)

g1 <- grad(obj, phi, y=y, tau=tau, D=D, lambda=lambda)
g2 <- grad_obj(phi, y, tau, D, lambda)

H1 <- genD(grad_obj, phi, y=y, 
              tau=tau, D=D, lambda=lambda, 
              )
H2 <- hess_obj(phi, y, tau, D, lambda)

tmp <- seq(-1, 1, .01)
h1 <- approx_checkLoss(tmp, .1)

f1 <- sapply(tmp, check_loss, tau=tau)
f2 <- approx_checkLoss(tmp, tau)
plot(f2~tmp, type="l")
lines(f1~tmp, col="red")
f2[201]
f1[201]
###############################################################################

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


