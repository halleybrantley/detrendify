library(devtools)
library(splines)
library(ggplot2)
library(microbenchmark)
library(kernlab)
load_all("detrendr")
rm(list=ls())
source("sim_generating_functions.R")
set.seed(12345678)
tau <- c(0.05, .1)
simDesign <- "peaks"
df <- generate_peaks_design(100000)
n <- 20000
overlap <- 2000
y <- df$y[1:n]

lambda = n
k=3
window_size = as.integer(n/2+overlap/2)
rho <- 1
max_iter = 4
quad = TRUE
use_gurobi = TRUE

trend0 <- get_trend(y, tau, lambda, k)

trend <- get_trend_windows(y, tau, lambda, k=3, 
                  window_size, 
                  overlap = overlap, max_iter = 4, quad = TRUE, 
                  use_gurobi = TRUE, rho = 1, update = 1)

plot(y, type="l")
lines(trend[,1], col="red")
lines(trend0[,1], col="blue")
max(abs(trend - trend0))
sd(y)
