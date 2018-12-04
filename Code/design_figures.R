library(ggplot2)
library(tidyverse)
source("trueQuantile.R")

load("../SimData/shapebeta_n_500_sim027.RData")
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
for (i in 1:length(tau)){
  df[,paste0("q",i)] <- trueQuantile("shapebeta", df$x, tau[i])
}
df_long <- df %>% gather("quantile", "trend", -c(y,x,f))

ggplot(df, aes(x=x, y=y)) + 
  geom_point(col="grey") +
  theme_bw() +
  geom_line(data=df_long, aes(x=x, y=trend, col=quantile)) +
  scale_color_brewer(palette = "Set1", 
                     labels = tau) +
  labs(title = "Beta")
ggsave("../Manuscript/Figures/shapebeta.png", width = 3, height = 3)  

load("../SimData/gaus_n_500_sim027.RData")
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
for (i in 1:length(tau)){
  df[,paste0("q",i)] <- trueQuantile("gaus", df$x, tau[i])
}
df_long <- df %>% gather("quantile", "trend", -c(y,x,f))

ggplot(df, aes(x=x, y=y)) + 
  geom_point(col="grey") +
  theme_bw() +
  geom_line(data=df_long, aes(x=x, y=trend, col=quantile)) +
  scale_color_brewer(palette = "Set1", 
                     labels = tau) +
  labs(title = "Gaussian")
ggsave("../Manuscript/Figures/gaus.png", width = 3, height = 3) 

load("../SimData/mixednorm_n_500_sim027.RData")
for (i in 1:length(tau)){
  df[,paste0("q",i)] <- trueQuantile("mixednorm", df$x, tau[i])
}
df_long <- df %>% gather("quantile", "trend", -c(y,x,f))

ggplot(df, aes(x=x, y=y)) + 
  geom_point(col="grey") +
  theme_bw() +
  geom_line(data=df_long, aes(x=x, y=trend, col=quantile)) +
  scale_color_brewer(palette = "Set1", 
                     labels = tau) +
  labs(title = "Mixed Normal")
ggsave("../Manuscript/Figures/mixednorm.png", width = 3, height = 3) 
