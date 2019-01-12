rm(list=ls())
load("../TimingData/window_times.RData")
all.times0 <- all.times
# load("../TimingData/window_times_2.RData")
# all.times <- rbind(all.times, all.times0)
windows.df <- data.frame(n = data_lengths,
                      median = apply(all.times, 2, median)*1e-9,
                      mean = apply(all.times, 2, mean)*1e-9)

load("../TimingData/single_window_times.RData")
all.times0 <- all.times
load("../TimingData/single_window_times_2.RData")
all.times0 <- rbind(all.times, all.times0)
# load("../TimingData/single_window_times_3.RData")
# all.times <- rbind(all.times, all.times0)

data_lengths <- seq(2000, 50000, 2000)
time.df <- data.frame(n = data_lengths, 
                      median = apply(all.times, 2, median, na.rm=T)*1e-9,
                      mean = apply(all.times, 2, mean, na.rm=T)*1e-9)

ggplot(time.df, aes(x=n, y=median)) + 
  geom_line() + 
  geom_point() +
  geom_line(data = windows.df, col="red") +
  geom_point(data = windows.df, col="red") +
  theme_bw() +
  labs(y = "Time (s)", x = "Data Size")
ggsave("../Manuscript/Figures/Fig_timing_experiment.png", width = 7, height = 2)
