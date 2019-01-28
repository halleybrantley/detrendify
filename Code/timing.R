library(devtools)
library(splines)
library(ggplot2)
library(microbenchmark)
load_all("detrendr")
rm(list=ls())
source("sim_generating_functions.R")
tau <- c(0.05, 0.1, 0.2)
simDesign <- "peaks"


times <- 25
overlap <- 500
data_lengths <- seq(6000, 51000, 3000)
times_df <- data.frame(n = NA, t = NA, time = NA, windows = NA)
diff_df <- data.frame(n = NA, t = NA, maxdiff = NA, sd = NA)
i <- 1
j <- 1
# load("../TimingData/timing_comparison.RData")
# i <- nrow(times_df)
# j <- nrow(diff_df)
for (n in data_lengths){
  for (t in 1:times){
    print(sprintf("n = %i", n))
    print(sprintf("t = %i", t))
    df <- generate_peaks_design(n)
    for (w in 1:4){
      times_df$n[i] <- n
      times_df$t[i] <- t
      times_df$windows[i] <- w

      if (w == 1) {
        times_df$time[i] <- microbenchmark(
                      trend <- get_trend(df$y, tau, lambda = n/5, k=3), 
                      times = 1)$time*1e-9 
      } else {
        times_df$time[i] <- 
          microbenchmark( 
            trend_w <- get_trend_windows(df$y, tau, lambda = n/5, k=3, 
                                 window_size = round(n/w+overlap*(w-1)/w), 
                                 overlap = overlap, max_iter = 20, 
                                 update = 10, rho = 3, eps_abs = .01, 
                                 eps_rel = 1e-3, 
                                 quad = TRUE, use_gurobi = TRUE), 
            times = 1)$time*1e-9
        diff_df$n[j] <- n
        diff_df$t[j] <- t
        diff_df$windows[j] <- w
        diff_df$maxdiff[j] <- max(abs(trend_w - trend))
        diff_df$sd[j] <- sd(df$y)
        j <- j+1
        diff_df[j,] <- NA
      }
      i <- i+1
      times_df[i,] <- NA
	save(times_df, diff_df, file = "../TimingData/timing_comparison.RData")
    }
  }
}

# trend <- get_trend_windows(df$y, tau, lambda = n/5, k=3, 
#                   window_size = round(n/w+overlap*(w-1)/w), 
#                   overlap = overlap, max_iter = 20, 
#                   update = 1, rho = 3, eps_abs = .01, 
#                   eps_rel = 1e-3, 
#                   quad = TRUE, use_gurobi = TRUE)
# plot(trend[,3], type='l')

