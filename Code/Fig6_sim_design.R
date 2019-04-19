################################################################################
# Simulation Design Example Figures
################################################################################
library(ggplot2)
library(tidyverse)
rm(list=ls())
source("sim_generating_functions.R")

load("../SimData/shapebeta_n_500_sim027.RData")
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
for (i in 1:length(tau)){
  df[,paste0("q",i)] <- trueQuantile("shapebeta", df$x, tau[i])
}
df_long <- df %>% gather("quantile", "trend", -c(y,x,f))
text_size <- 20

ggplot(df, aes(x=x, y=y)) + 
  geom_point(col="grey") +
  theme_bw() +
  geom_line(data=df_long, aes(x=x, y=trend, col=quantile)) +
  scale_color_brewer(palette = "Set1", 
                     labels = tau) +
  labs(title = "Beta", x = "", y="") +
  guides(col = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-.05, 1.05))+
  theme(text = element_text(size=text_size), 
        plot.title = element_text(size = text_size))
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
  labs(title = "Gaussian", x = "", y="") +
  guides(col = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-.05, 1.05))+
  theme(text = element_text(size=text_size), 
        plot.title = element_text(size = text_size))
ggsave("../Manuscript/Figures/gaus.png", width = 2.8, height = 3) 

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
  labs(title = "Mixed Normal", x = "", y="", col = "Quantile") +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(-.05, 1.05))+
  theme(text = element_text(size=text_size), 
        plot.title = element_text(size = text_size), 
        legend.title = element_text(size = (text_size-2)))
ggsave("../Manuscript/Figures/mixednorm.png", width = 4.5, height = 3) 

n <- 1000
i <- 5
simDesign <- "peaks"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
ggplot(df, aes(x=x, y=y, col = "data")) + 
  geom_line() +
  geom_line(aes(y=peaks+baseline, col = "signal + baseline")) +
  geom_line(aes(y=baseline, col = "baseline")) +
  theme_bw() + 
  scale_color_manual(values = c("red", "grey", "blue")) + 
  labs(x = "t", col = "")
ggsave("../Manuscript/Figures/ex_peaks.png", width = 7, height = 2.5)
