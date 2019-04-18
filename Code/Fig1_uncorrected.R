################################################################################
# Figure 1 - full day of uncorercted data on 4-13-17
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
rm(list=ls())
colPal <- c('#1b7837', '#762a83')
tau <- c(0.01, 0.05, 0.1)
nodes <- c("c", "d", "e")
load(sprintf("../SPod/spodPIDs.RData"))
spod_trends <- data.frame(time = spodPIDs$time)

for (node in nodes){
  load(file=sprintf("../SPod/node_%s_trend.RData", node))
  spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
  names(spod_trends)[(ncol(spod_trends)-length(tau)+1):ncol(spod_trends)] <-
    paste(node, tau, sep = "_")
}

spodRaw <- spodPIDs %>%
  gather("node", "PID", -time)

ggplot(spodRaw, aes(x=time, y=PID, col=node)) +
  geom_line(alpha=0.5) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", labels = c("a", "b", "c"))+
  labs(col="SPod")

ggsave("../Manuscript/Figures/uncorrected_data.png", width = 7, height = 2.5)
