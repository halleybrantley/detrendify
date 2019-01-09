library(devtools)
library(splines)
library(ggplot2)
library(microbenchmark)
library(kernlab)
load_all("detrendr")
rm(list=ls())
source("sim_generating_functions.R")
set.seed(12345678)
tau <- c(0.05)
simDesign <- "peaks"
df <- generate_peaks_design(1000)
n <- 20
overlap <- 5
y <- df$y[1:n]

trend <- get_trend_windows(y, tau, lambda = n, k=3, 
                  window_size = n/2+overlap/2, 
                  overlap = overlap, max_iter = 4, quad = TRUE, 
                  use_gurobi = FALSE, rho = .5)
plot(y, type="l")
lines(trend[,2], col="blue")
