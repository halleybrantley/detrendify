################################################################################
# Peaks Simulation Results
# Halley Brantley
################################################################################
library(tidyverse)
rm(list=ls())
i <- 1
n <- 500
tau <- c(0.01, 0.05, 0.1)
simDesign <- "peaks"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData",
             "detrend_eBIC", simDesign, n, i))
trend_detrend <- as.data.frame(trend)
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData",
             "windows", simDesign, n, i))
trend_qsreg <- as.data.frame(trend)

trend_detrend$x <- df$x
trend_detrend$q <- df$baseline+qnorm(tau[3], sd = 0.25)
trend_detrend$q2 <- df$baseline+qnorm(tau[2], sd = 0.25)
trend_detrend$qsreg <- trend_qsreg[,3]
trend_detrend$qsreg2 <- trend_qsreg[,2]

ggplot(df, aes(x=x, y=y)) +
  geom_line(col="grey") +
  geom_line(data=trend_detrend, aes(y=q, x=x, col = "true 0.1"), lwd = 1.5) +
  geom_line(data=trend_detrend, aes(y=q2, x=x, col = "true 0.05"), lwd = 1.5) +
  #geom_line(data=trend_detrend, aes(y=V3, x=x, col = "detrendr 0.1")) +
  geom_line(data=trend_detrend, aes(y=qsreg, x=x, col = "qsreg 0.1")) +
  geom_line(data=trend_detrend, aes(y=V2, x=x, col = "detrendr 0.05")) +
  geom_line(data=trend_detrend, aes(y=qsreg2, x=x, col = "qsreg 0.05")) +
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(col="")
ggsave("../Manuscript/Figures/ex_baseline.png", width = 7, height = 4)
