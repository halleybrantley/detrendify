library(devtools)
library(splines)
library(ggplot2)
library(microbenchmark)
load_all("detrendr")
rm(list=ls())
source("sim_generating_functions.R")
tau <- c(0.05, 0.1)
simDesign <- "peaks"


plot(y~x, df, type="l")
lines((peaks+baseline)~x, df, col="blue")
lines((baseline)~x, df, col="red")


times <- 3
overlap <- 1000
data_lengths <- seq(50000, 65000, 5000)
single.times <- matrix(NA, ncol = length(data_lengths), nrow = times)
window2.times <- matrix(NA, ncol = length(data_lengths), nrow = times)
window3.times <- matrix(NA, ncol = length(data_lengths), nrow = times)
i <- 1
for (n in data_lengths){
  print(sprintf("n = %i", n))
  df <- generate_peaks_design(n)
  single.times[,i] <- microbenchmark(
    trend <- get_trend(df$y, tau, lambda = c(n,n), k=3), 
    times = times)$time
  window2.times[,i] <- microbenchmark(
    trend_w <- get_trend_windows(df$y[1:n], tau, lambda = c(n,n), k=3, 
                                 window_size = as.integer(n/2+overlap/2), 
                                 overlap = overlap, max_iter = 4, quad = TRUE, 
                                 use_gurobi = TRUE), 
    times = times)$time
  window3.times[,i] <- microbenchmark(
    trend_w3 <- get_trend_windows(df$y[1:n], tau, lambda = c(n,n), k=3, 
                                 window_size = round(n/3+2*overlap/3) + 1, 
                                 overlap = overlap, max_iter = 4, quad = TRUE, 
                                 use_gurobi = TRUE), 
    times = times)$time
  i <- i+1
  
  save(df, trend, trend_w, trend_w3, 
       file = sprintf("../TimingData/trends_%i.RData", n))
  save(single.times, window2.times, window3.times, data_lengths, 
       file="../TimingData/all_times.RData")
}

# load("../TimingData/single_window_times_3.RData")
# time.df <- data.frame(n = data_lengths, 
#            median = apply(all.times, 2, median)*1e-9,
#            mean = apply(all.times, 2, mean)*1e-9)
# 
# ggplot(time.df, aes(x=n, y=median)) + geom_line() + geom_point()
# 

