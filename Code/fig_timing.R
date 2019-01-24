library(tidyverse)
rm(list=ls())

load("../TimingData/timing_comparison.RData")
times_df <- times_df[-nrow(times_df),]    
diff_df <- diff_df[-nrow(diff_df),]  


times_df %>% group_by(n, windows) %>%
  summarise(time_mean = mean(time), 
            se_time = sd(time)/sqrt(n())) %>%
  ggplot(aes(x=n, y=time_mean, group=factor(windows), col=factor(windows))) + 
  geom_line() + 
  geom_linerange(aes(ymin = time_mean - 2*se_time, ymax = time_mean + 2*se_time))+
  theme_bw() +
  scale_color_brewer(palette="Set1")+
  labs(y = "Time (s)", col = "# Windows")
ggsave("../Manuscript/Figures/Fig_timing_experiment.png", width = 6, height = 2.5)

diff_df %>% group_by(windows) %>%
  summarise(max = max(maxdiff), 
            mean = mean(maxdiff),
            median = median(maxdiff))
