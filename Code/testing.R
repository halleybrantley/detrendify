require(devtools)
require(Matrix)
require(ggplot2)
library(lpSolve)
library(rbenchmark)
library(Rglpk)
library(microbenchmark)
library(fasta)

rm(list=ls())
load_all("../Code/detrendr")

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
spodNode <- spodNode[35001:38600,]
plot(pid~time, spodNode, type="l", main = "raw")


y <- spodNode$pid/1000
k <- 3
tau <- .05
D <- get_Dk(length(y), k)
lambda <- 5

# benchmark(replications = 3, 
#           glpk = lpglpk_trendfilter(y, tau, lambda, D),
#           cvx = cvx_trendfilter(y, tau, lambda, D), 
#           columns=c('test', 'replications', 'elapsed'))


theta_lpgl <- lpglpk_trendfilter(y, tau, lambda, D)
plot(pid/1000~time, spodNode, type="l", main = "raw")
lines(theta_lpgl~spodNode$time, col="red")

eval_obj <- function(theta, y, tau, eta, lambda) {
  val <- sum((0.5 * abs(y-theta) + (tau - 0.5) * (y-theta))/length(y),
             lambda*norm(eta,"o"))
  return(val)
}
eval_obj(theta_lpgl, y, tau, D%*%theta_lpgl, lambda)

eta_lpgl <- D%*%theta_lpgl

x <- c(theta_lpgl, as.numeric(eta_lpgl))

nu <- .01
n <- length(y)
step_size <- 1e-4
f(x, y, tau, D, nu) 
gf <- gradf(x, y, tau, D, nu)
pg <- proxg(x, D, lambda, step_size)
plot(pg[1:n]~x[1:n])
plot(pg[(n+1):length(x)]~x[(n+1):length(x)])

x_f <- fasta(f, gradf, g, proxg)
  
x <- x - gf*.01
f(x, y, tau, D, nu) 

m <- length(nrow(D))
theta <- y
eta <- as.numeric(D%*%theta)


M <- Cholesky(Diagonal(n) + Matrix::crossprod(D))
f0 <- eval_obj(theta, y, tau, D%*%theta, lambda)
iters <- 10000
fvals <- double(iters)
gaps <- double(iters)
step <- double(iters)


# Spingarn

for (i in 1:iters){
  theta_old <- theta
  eta_old <- eta
  prox_sol <- prox_R(theta, eta, y, lambda, tau, 1)
  theta <- prox_sol$theta
  eta <- prox_sol$eta
  proj_sol <- project_V_R(2*theta - theta_old, 2*eta - eta_old, D, M)
  theta <- theta_old + 1*(proj_sol$theta - theta)
  eta <- eta_old + 1*(proj_sol$eta - eta)
  f <- eval_obj(proj_sol$theta, y, tau, proj_sol$eta, lambda)
  fvals[i] <- f
}

theta_s <- prox_R(theta, eta, y, lambda, tau, 1)$theta

n <- length(y)
m <- length(nrow(D))
theta <- y
eta <- as.numeric(D%*%theta)
mu <- 100
tau1 <- 2
tau2 <- 1.9
primals <- double(iters)
duals <- double(iters)
rho <- 1
u_theta <- double(n)
u_eta <- double(m)

# Scaled ADMM Algorithm
for(i in 1:iters){
  prox_sol <- prox_R(proj_sol$theta - u_theta, proj_sol$eta - u_eta, y, lambda, 
                     tau, 1/rho)
  proj_sol <- project_V_R(prox_sol$theta + u_theta, prox_sol$eta + u_eta, D, M)

  u_theta <- u_theta + prox_sol$theta - proj_sol$theta
  u_eta <- u_eta + prox_sol$eta - proj_sol$eta
  
  primal_resid <- norm(prox_sol$eta - D%*%prox_sol$theta, "f")
  #/(max(norm(prox_sol$eta, "f"), norm(D%*%prox_sol$theta, "f")))
  dual_resid <- norm(rho*crossprod(D, eta-proj_sol$eta), "f") 
  #/ norm(crossprod(D, u_eta), "f")
  
  eta <- proj_sol$eta
  
  f <- eval_obj(proj_sol$theta, y, tau, proj_sol$eta, lambda)
  

  if (primal_resid > mu*dual_resid){
    rho <- rho*tau1
    u_theta <- u_theta/tau1
    u_eta <- u_eta/tau1
    print(rho)
  } else if (dual_resid > mu*primal_resid) {
    rho <- rho/tau2
    u_theta <- u_theta*tau2
    u_eta <- u_eta*tau2
    print(rho)
  }

  primals[i] <- primal_resid
  duals[i] <- dual_resid
  fvals[i] <- f

}

plot(y, type="l")
lines(theta_s, type="l", col="red")
lines(prox_sol$theta, type="l", col="blue")
lines(theta_lpgl, col="green")

lines(theta_cvx, col="blue")

plot(primals[4000:iters], type="l")
plot(duals[4000:iters], type="l")

resids <- primals + duals
plot(resids[4000:5000], type="l")

plot(duals[100:500], type="l")

plot(fvals[4000:iters], type="l")
min(fvals)
sol$obj


###############################################################################

n <- length(y)
m <- length(nrow(D))
theta <- y
M <- Cholesky(Diagonal(n) + Matrix::crossprod(D))
f0 <- eval_obj(theta, y, tau, D%*%theta, lambda)
iters <- 5500
fvals <- double(iters)
gaps <- double(iters)
step <- double(iters)
beta <- .995


L <- 1
mu <- L/norm(D, "f")^2
z <- prox_f2_R(as.numeric(D%*%theta), lambda, L)
u <- D%*%theta - z

for(i in 1:iters){
  theta <- prox_f1_R(as.numeric(theta - mu/L * crossprod(D, D %*% theta - z + u)), 
                     y, tau, mu)
  z <- prox_f2_R(as.numeric(D%*%theta + u), lambda, L)
  u <- u + D%*%theta - z
  f <- eval_obj(theta, y, tau, D %*% theta, lambda)
  gaps[i] <- gap
  fvals[i] <- f
  step[i] <- steps[1]
}

plot(y, type="l")
lines(theta, type="l", col="blue")
lines(theta_cvx, col="red")
plot(fvals, type="l")
