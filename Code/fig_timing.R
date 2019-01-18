library(tidyverse)
rm(list=ls())
load("../TimingData/stopping_crit_times.RData")
all.times0 <- as.data.frame(t(all.times))
all.times0$n <- data_lengths

load("../TimingData/stopping_crit_times_2.RData")
all.times1 <-  as.data.frame(t(all.times))
all.times1$n <- data_lengths
windows.df <- bind_rows(all.times0, all.times1) %>% 
  gather("rep", "time", -n) %>%
  group_by(n) %>%
  summarise(time = mean(time, na.rm=T)*1e-9)


load("../TimingData/single_window_times_4.RData")
all.times0 <- as.data.frame(t(all.times))
all.times0$n <- data_lengths

load("../TimingData/single_window_times.RData")
all.times <- as.data.frame(t(all.times))
all.times$n <- data_lengths


time.df <- bind_rows(all.times, all.times0) %>% 
  gather("rep", "time", -n) %>%
  group_by(n) %>%
  summarise(time = median(time, na.rm=T)*1e-9)

ggplot(time.df, aes(x=n, y=time)) + 
  geom_point() +
  geom_line() +
  geom_line(data = windows.df, col="red", aes(y=time)) +
  geom_point(data = windows.df, col="red", aes(y=time)) +
  theme_bw() +
  labs(y = "Time (s)", x = "Data Size")
ggsave("../Manuscript/Figures/Fig_timing_experiment.png", width = 7, height = 2)
