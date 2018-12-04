library(devtools)
library(rbenchmark)
library(quantreg)
library(fields)
library(tidyverse)
load_all("detrendr")
rm(list=ls())
source("data_generating_functions.R")
source("npqw/lib_npqw_mixed.R")


tau <- .1
n <- 500
overlap <- 300
timing <- {}
for (n in seq(500, 5000, 500)){
  sim_data <- generate_gaus(n)
  overlap <- n/10
  time_n <- benchmark(
    "single" = {
      f_detrend <- gurobi_trend(sim_data$y, tau, lambda=n^2, k=3)
    },
    "consensus" = {
      f_consensus <- consensus_ADMM(y=sim_data$y, tau, lambda=n^2, k=3, rho=10, 
                                    window_size=(n+overlap)/2, overlap=overlap, 
                                    max_iter=2, eps = 2e-5, update = 100)$theta
    },  "qsreg" = {
      f_qsreg <- predict(qsreg(sim_data$x, sim_data$y, lam=1e-11, alpha=tau))
    },
    # "npqw" = {
    #   sd.pilot <- sd.fan.yao(x=sim_data$x,
    #                          y=sim_data$y,
    #                          regtype=regtype,
    #                          ckertype=ckertype,
    #                          nmulti=nmulti)
    #   Q <- qnorm(0.1,mean=sim_data$y,sd=sd.pilot)
    #   f_npqw <- fitted(npreg(txdat=sim_data$x,
    #                          tydat=Q,
    #                          bws=0.01,
    #                          regtype=regtype,
    #                          okertype="liracine",
    #                          ukertype="liracine",               
    #                          ckertype=ckertype))
    # }, 
    "rqss" = {
      f_rqss <- predict(rqss(y ~ qss(x, lambda = .2), tau = tau, 
                             data = sim_data), sim_data)
    },
    replications = 1)

  print(paste("n =", n))
  time_n$n <- n
  timing <- rbind(timing, time_n[, c("n", "test", "elapsed")])
}


ggplot(timing, aes(x=n, y=elapsed, col=test)) + geom_line() 

regtype <- "ll"
ckertype <- "epanechnikov"
nmulti <- 1

timing2 <- {}
for (n in seq(500, 10000, 500)){
  sim_data <- generate_gaus(n)
  # Choose qsreg tuning parameter
  # fit.cv <- qsreg(sim_data$x, sim_data$y, alpha=tau)
  # lam <- fit.cv$cv.grid[which.min(fit.cv$cv.grid[,"CV"]), "lambda"]
  # 
  # Chose npqw tuning parameter

  
  # optim.out <- npqw.optim(x=sim_data$x,
  #                         y=sim_data$y,
  #                         tau=tau,
  #                         nmulti=nmulti,
  #                         sd.pilot=sd.pilot,
  #                         ckertype=ckertype,
  #                         regtype=regtype)
  # 
  #   delta <- optim.out$delta
  #   bws <- optim.out$bws
    
  time_n <- benchmark(
    "qsreg" = {
      f_qsreg <- predict(qsreg(sim_data$x, sim_data$y, lam=1e-11, alpha=tau))
    },
    "npqw" = {
      sd.pilot <- sd.fan.yao(x=sim_data$x,
                             y=sim_data$y,
                             regtype=regtype,
                             ckertype=ckertype,
                             nmulti=nmulti)
      Q <- qnorm(0.1,mean=sim_data$y,sd=sd.pilot)
      f_npqw <- fitted(npreg(txdat=sim_data$x,
                      tydat=Q,
                      bws=0.01,
                      regtype=regtype,
                      okertype="liracine",
                      ukertype="liracine",               
                      ckertype=ckertype))
    }, 
    "rqss" = {
      f_rqss <- predict(rqss(y ~ qss(x, lambda = .2), tau = tau, 
                   data = sim_data), sim_data)
    },
    replications = 3)
  
  print(paste("n =", n))
  timing2 <- rbind(timing2, c(n, time_n$elapsed))
}
