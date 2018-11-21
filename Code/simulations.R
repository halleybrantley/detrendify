# Simulations
library(fields)
library(quantreg)
library(devtools)
# library(qrsvm)
library(np)
library(microbenchmark)
load_all("detrendr")
library(tidyverse)
#source("function_jcgs.R")
rm(list=ls())
nSim <- 100
n <- 500
x <- (seq(1,n,1)-1)
tau <- .1
set.seed(986532)
f_true <- sin(10*x/n) + (x/n + .25)*qnorm(tau, 0, .1)/.1
loss_df <- data.frame(SimNum = 1:nSim, 
                      trend_SIC = NA, trend_valid = NA, rqss = NA, 
                      qsreg = NA, npqreg = NA)

for (i in 1:nSim){
  y <- sin(10*x/n) + (x/n + .25)*rnorm(n, 0, .1)/.1
  lam_SIC <- lambda_SIC(y, tau, 3)
  lam_valid <- lambda_valid(y, tau, 3, 5)
  lam_SIC2 <- lambda_SIC(y, tau, 2, 
                         lambdaSeq = seq(.5, length(y)/10, .2))
  
  f_trend_SIC <- gurobi_trend(y, tau, lam_SIC$lambda, k=3)
  loss_df$trend_SIC[i] <- mean((f_trend_SIC-f_true)^2)
  
  f_trend_valid <- gurobi_trend(y, tau, lam_valid$lambda, k=3)
  loss_df$trend_valid[i] <- mean((f_trend_valid-f_true)^2)
  
  fit_rqss <- rqss(y ~ qss(x, lambda = 2*lam_SIC2$lambda), tau = tau)
  f_rqss <- predict(fit_rqss, data.frame(x=x))
  loss_df$rqss[i] <- mean((f_rqss-f_true)^2)
  
  fit_qsreg <- qsreg(x, y, maxit.cv = 50, alpha=tau, hmin = -9)
  f_qsreg <- predict(fit_qsreg)
  loss_df$qsreg[i] <- mean((f_qsreg-f_true)^2)
  
  bw <- npcdistbw(y~x)
  sink("NUL")
  fit_npqreg <- npqreg(bws=bw, tau=tau)
  sink()
  f_npqreg <- fitted(fit_npqreg)
  loss_df$npqreg[i] <- mean((f_npqreg-f_true)^2)
  print(sprintf("\nSim %i complete", i))
}
  
plot(f_true~x, type="l")
lines(f_rqss~x, col="red")
lines(f_trend_SIC~x, col="blue")
lines(f_trend_valid~x, col="purple")
lines(f_qsreg~x, col="darkgreen")
lines(f_npqreg~x, col="cyan")

loss_long <- gather(na.omit(loss_df), "method", "mse", -SimNum)
ggplot(loss_long, aes(x=method, y=mse)) + geom_boxplot() +
  stat_summary(fun.data = "mean_cl_boot", col="red") +
  geom_hline(yintercept = mean(loss_df$trend_SIC))

timing <- {}
for (n in seq(100, 2000, 200)){
  x <- (seq(1,n,1)-1)
  y <- sin(10*x/n) + (x/n + .25)*rnorm(n, 0, .1)/.1
  lam_SIC <- lambda_SIC(y, tau, 3)
  lam_SIC2 <- lambda_SIC(y, tau, 2, 
                         lambdaSeq = seq(.5, length(y)/10, .2))
  #bw <- npcdistbw(y~x)
  
  time_n <- summary( microbenchmark(
    f_trend_valid <- gurobi_trend(y, tau, lam_valid$lambda, k=3),
    fit_rqss <- rqss(y ~ qss(x, lambda = 2*lam_SIC2$lambda), tau = tau),
    fit_qsreg <- qsreg(x, y, maxit.cv = 50, alpha=tau, hmin = -9), 
    # {sink("NUL")
    #  fit_npqreg <- npqreg(bws=bw, tau=tau)
    #  sink()
    # }, 
    times = 20)) %>% 
      select(expr, mean)
  time_n$n <- n
  timing <- rbind(timing, time_n)
}

ggplot(timing, aes(x=n, y=mean, col=expr)) + geom_line() +
  ylim(c(0,500))

mean((f_true-f_trend)^2)
mean((f_true - f_trend_SIC)^2)
mean((f_true - f_trend_valid)^2)
mean((f_true - f_qsreg)^2)



f <- predict(fit_qsreg)
resid <- y-f
w <- rep(tau, length(resid))
w[resid < 0] <- tau-1
w <- w/(2*resid)
w <- w/sum(w)

fit_sq <- sreg(x, y, weights = w)
fit_sq2 <- smooth.spline(x,y, w)
A <- fit_sq$diagA
f2 <- predict(fit_sq)
resid <- y-f2

LOOCV_mse <- mean(resid^2/(1-A)^2)
# GCV MSE with weights
GCV_mse <- crossprod(resid, diag(w))%*%resid/(1-mean(A))^2/n

summary(fit_sq)


# Fit Koenker Smoothing Splines
rmse_qss <- function(x, y, tau, lambda, f_true){
  fit <- rqss(y ~ qss(x, lambda = lambda), tau = tau)
  fhat <- fit$coef[1] + fit$coef[-(1:2)]
  rmse <- sqrt(mean((f_true[-(1:2)] - fhat)^2))
  return(rmse)
}


# Fit trendfiltering
rmse_trendFilter <- function(y, x, tau, k, lambda, f_true) {
  f_trend <- gurobi_trend(y, tau, lambda, k=k)
  rmse <- sqrt(mean((f_true[-(1:2)] - f_trend[-(1:2)])^2))
  return(rmse)
}

# Fit Method by Oh and Nychka
rmse_Yu <- function(y, x, tau, f_true){
  fit_Yu <- aresult(x,(2000-100+1),0.04,hry1(100,2000,x),hry2(100,tau,2000,y)) #Yu's method
  rmse <- sqrt(mean((f_true[-(1:2)] - fit_Yu[-(1:2)])^2))
  return(rmse)
}

# Fit Method by Oh and Nychka
rmse_ES <- function(y, x, tau, f_true){
  finit <- aresult(x,(2000-100+1),0.04,hry1(100,2000,x),hry2(100,tau,2000,y)) #Yu's method
  fit_Oh <- qreq1d.sreg(x,y,finit,tau,0.2) # the proposed method
  rmse <- sqrt(mean((f_true[-(1:2)] - fit_Oh[-(1:2)])^2))
  return(rmse)
}

# Generate Data
y <- sin(10*x) + (x + .25)*rnorm(2000, 0, .07)/.1
# y2 <- sin(10*x) + rgamma(2000, 3, 1)

# True Quantile
tau <- 0.1
f_true <- sin(10*x) + (x + .25)*qnorm(tau, 0, .07)/.1


# Choose Lamda
sic <- function(lam, x, y, tau){
  # -1 indicated Schwarz criteria k = log(n)
  a <- AIC(rqss(y ~ qss(x, lambda = lam), tau=tau, method = "sfn"), k=-1)
}


lam_qss <- optimize(sic, c(.0001, 10), x=x, y=y, tau = tau)
folds <- sample(1:4, n, replace = TRUE)
lam_trend_k3 <- optimize(valid_checkLoss, c(n/3, 3*n), folds=folds, y=y, tau=0.5, k=3)
lam_trend_k2 <- optimize(valid_checkLoss, c(n/3, 3*n), folds=folds, y=y, tau=0.5, k=2)

set.seed(987654)

rmse_comp <- data.frame(i = 1:50, QSS = NA, ES = NA, Yu = NA, trend3 = NA)

for (i in 1:50){
  y <- sin(10*x) + (x + .25)*rnorm(2000, 0, .07)/.1
  rmse_comp[i,"ES"] <- rmse_ES(y, x, tau, f_true)
  rmse_comp[i, "QSS"] <- rmse_qss(x, y, tau, lam_qss$minimum, f_true)
  rmse_comp[i, "Yu"] <- rmse_Yu(y, x, tau, f_true) 
  rmse_comp[i, "trend3"] <- rmse_trendFilter(y, x, tau, k=3, lam_trend_k3$minimum, f_true) 
}

rmse_long <- rmse_comp %>% gather("method", "rmse", -i)
colMeans(rmse_comp)
apply(rmse_comp, 2, sd)/sqrt(50)

ggplot(rmse_long, aes(x=method, y=rmse)) + geom_boxplot()

plot(y~x, col="grey")
lines(x[-(1:2)],fhat, col="blue")
lines(f_true~x)
lines(f_trend~x, col="red")
lines(fit_Oh~x, col="purple")

set.seed(112218)
y <- sin(10*x) + rgamma(2000, 3, 1)
lam_qss <- optimize(sic, c(.0001, 10), x=x, y=y, tau = tau)
folds <- sample(1:5, n, replace = TRUE)
lam_trend_k3 <- optimize(valid_checkLoss, c(n/3, 3*n), folds=folds, y=y, tau=0.5, k=3)
lam_trend_k2 <- optimize(valid_checkLoss, c(n/3, 3*n), folds=folds, y=y, tau=0.5, k=2)

tau <- .9
rmse_gamma <- data.frame(i = 1:50, QSS = NA, ES = NA, trend2 = NA, trend3 = NA)
f_true <- sin(10*x) + qgamma(tau, 3, 1)

for (i in 1:50){
  y <- sin(10*x) + rgamma(2000, 3, 1)
  rmse_gamma[i,"ES"] <- rmse_ES(y, x, tau, f_true)
  rmse_gamma[i, "QSS"] <- rmse_qss(x, y, tau, lam_qss$minimum, f_true)
  rmse_gamma[i, "trend2"] <- rmse_trendFilter(y, x, tau, k=2, lam_trend_k2$minimum, f_true) 
  rmse_gamma[i, "trend3"] <- rmse_trendFilter(y, x, tau, k=3, lam_trend_k3$minimum, f_true) 
}

rmse_gamma_long <- rmse_gamma %>% gather("method", "rmse", -i)
colMeans(rmse_gamma)
apply(rmse_gamma, 2, sd)/sqrt(50)

ggplot(rmse_gamma_long, aes(x=method, y=rmse)) + geom_boxplot()

plot(y~x, col="grey")
lines(x[-(1:2)],fhat, col="blue")
lines(f_true~x)
lines(f_trend~x, col="red")
lines(fit_Oh~x, col="purple")
