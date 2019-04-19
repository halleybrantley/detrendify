################################################################################
# Figure 9 - Peaks classification example
# Halley Brantley
################################################################################
library(ggplot2)
n <- 1000
i <- 5
method <- "detrend_eBIC"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", "peaks", n, i))
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             method, "peaks", n, i))
df$signal <- as.numeric(df$peaks > 0.5)
y_adj <- df$y - trend[,2]
signal_hat <- as.numeric(y_adj > 1.2)
trend <- as.data.frame(trend)
trend$x <- df$x
ggplot(df, aes(x=x, y=y)) + geom_line(col="grey") +
  geom_point(data=subset(df, signal==1), col="red", size = 1.4) +
  geom_point(data=df[signal_hat==1,], col="blue", size = 0.8) + 
  geom_line(data=trend, aes(y=V3)) +
  theme_bw() + 
  theme(text = element_text(size=16)) +
  labs(x="", y="")

ggsave("../Manuscript/Figures/peaks_eg_class.png", width = 7, height =2.5)
