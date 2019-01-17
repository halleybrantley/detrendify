library(tidyverse)
rm(list=ls())
load("../TimingData/stopping_crit_times.RData")
all.times0 <- all.times
# load("../TimingData/window_times_2.RData")
# all.times <- rbind(all.times, all.times0)
windows.df <- data.frame(n = data_lengths,
                      median = apply(all.times, 2, median)*1e-9,
                      mean = apply(all.times, 2, mean)*1e-9)


load("../TimingData/single_window_times_4.RData")
all.times <- as.data.frame(t(all.times))
all.times$n <- data_lengths
times.long <- all.times %>% gather("rep", "time", -n)
times.long$time <- times.long$time*1e-9

time.df <- data.frame(n = data_lengths, 
                      median = apply(all.times, 2, median, na.rm=T)*1e-9,
                      mean = apply(all.times, 2, mean, na.rm=T)*1e-9)

ggplot(time.df, aes(x=n, y=median)) + 
  geom_point() +
  geom_line() +
  geom_line(data = windows.df, col="red", aes(y=median)) +
  geom_point(data = windows.df, col="red", aes(y=median)) +
  theme_bw() +
  labs(y = "Time (s)", x = "Data Size")
ggsave("../Manuscript/Figures/Fig_timing_experiment.png", width = 7, height = 2)
