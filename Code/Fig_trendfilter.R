library(devtools)
install_github("hadley/l1tf")
install_github("glmgen/genlasso")
library(l1tf)
library(genlasso)
source("sim_generating_functions.R")
set.seed(1122334455)
n <- 500
df <- generate_mixednorm(n)
trend <- trendfilter(df$y, ord = 3)
df$trend <- coef(trend, lambda = 5000)$beta
df2 <- generate_peaks_design(n)
trend2 <- trendfilter(df2$y, ord = 3)
df2$trend <-  coef(trend2, lambda = 5000)$beta

ggplot(df, aes(x=x, y=y)) + 
  geom_line(col = "grey") + 
  geom_line(aes(y=trend), col="red") + 
  theme_bw()
ggsave("../Manuscript/Figures/trend_filter_eg1.png", width = 3, height = 2)

ggplot(df2, aes(x=x, y=y)) + 
  geom_line(col = "grey") + 
  geom_line(aes(y=trend), col="red") + 
  theme_bw()
ggsave("../Manuscript/Figures/trend_filter_eg2.png", width = 3, height = 2)

