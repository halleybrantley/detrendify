library(devtools)
library(splines)
library(ggplot2)
library(microbenchmark)
load_all("detrendr")
rm(list=ls())
source("sim_generating_functions.R")
set.seed(12345678)
tau <- c(0.05, 0.1)
simDesign <- "peaks"
df <- generate_peaks_design(100000)

plot(y~x, df, type="l")
lines((peaks+baseline)~x, df, col="blue")
lines((baseline)~x, df, col="red")


times <- 1
data_lengths <- seq(60000, 65000, 5000)
all.times <- matrix(NA, ncol = length(data_lengths), nrow = times)
trend_diff <- list()
overlap <- 2000
i <- 1
for (n in data_lengths){
  print(sprintf("n = %i", n))
  gc()
  all.times[,i] <- microbenchmark(
      trend_w <- get_trend_windows(df$y[1:n], tau, lambda = c(n,n), k=3, 
                                   window_size = as.integer(n/2+overlap/2), 
                                   overlap = overlap, max_iter = 4, quad = TRUE, 
                                   use_gurobi = TRUE), 
      times = times)$time
  trend <- get_trend(df$y[1:n], tau, lambda = c(n,n), k=3)
  
  trend_diff[[i]] <- trend_w - trend
  i <- i+1
  save(all.times, data_lengths, trend_diff,
       file="../TimingData/stopping_crit_times_2.RData")
}

# load("../TimingData/window_times.RData")
# time.df <- data.frame(n = data_lengths,
#            median = apply(all.times, 2, median)*1e-9,
#            mean = apply(all.times, 2, mean)*1e-9)
# 
# ggplot(time.df, aes(x=n, y=median)) + geom_line() + geom_point()
