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
df <- generate_peaks_design(50000)

plot(y~x, df, type="l")
lines((peaks+baseline)~x, df, col="blue")
lines((baseline)~x, df, col="red")


times <- 10
data_lengths <- seq(2000, 50000, 2000)
all.times <- matrix(NA, ncol = length(data_lengths), nrow = times)
i <- 1
for (n in data_lengths){
  print(sprintf("n = %i", n))
  gc()
  all.times[,i] <- microbenchmark(
    get_trend(df$y[1:n], tau, lambda = c(n,n), k=3), 
    times = times)$time
  i <- i+1
  save(all.times, data_lengths, 
       file="../TimingData/single_window_times_3.RData")
}

# load("../TimingData/single_window_times_3.RData")
# time.df <- data.frame(n = data_lengths, 
#            median = apply(all.times, 2, median)*1e-9,
#            mean = apply(all.times, 2, mean)*1e-9)
# 
# ggplot(time.df, aes(x=n, y=median)) + geom_line() + geom_point()
# 

