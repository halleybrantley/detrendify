rm(list=ls())
load("../TimingData/window_times.RData")
windows.df <- data.frame(n = data_lengths,
                      median = apply(all.times, 2, median)*1e-9,
                      mean = apply(all.times, 2, mean)*1e-9)

load("../TimingData/single_window_times.RData")
data_lengths <- seq(2000, 50000, 2000)
time.df <- data.frame(n = data_lengths, 
                      median = apply(all.times, 2, median)*1e-9,
                      mean = apply(all.times, 2, mean)*1e-9)
ggplot(time.df, aes(x=n, y=median)) + 
  geom_line() + 
  geom_point() +
  geom_line(data = windows.df, col="red") +
  geom_point(data = windows.df, col="red") 

