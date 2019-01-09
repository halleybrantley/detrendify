################################################################################
# Calculate mis-classification rates in peaks simulation
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
load_all("detrendr")
rm(list=ls())
tau <- c(0.01, 0.05, 0.1)
nSim <- 100
methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "qsreg", "rqss", "npqw") 
metrics <- data.frame(method = NA, n = NA, sim = NA, 
                      metric = NA, metric_type = NA, threshold = NA)

thresh <- c(.9, 1, 1.1, 1.2)
metric_types <- c("missclass", "precision", "recall", "CAA", "F1")
k <- 1

for (n in c(500,1000,2000,4000)){
  for (i in 1:nSim){
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", "peaks", n, i))
    df$signal <- as.numeric(df$peaks > 0.5)
    
    for (method in methods){
      load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
                   method, "peaks", n, i))
      for (threshold in thresh){
        for (metric in metric_types){
          for (j in 1:length(tau)){
            metrics$sim[k] <- i
            metrics$method[k] <- method
            metrics$n[k] <- n
            metrics$threshold <- threshold
            metrics$tau[k] <- tau[j]
            metrics$metric_type[k] <- metric
            metrics$metric[k] <- get_metric(trend[,j], df$y, df$signal, 
                                         threshold, metric)
            k <- k+1
            metrics[k, ] <- NA
          }
        }
      }
    }
  }
}
################################################################################
# Need to update

peaks_long <- miss_class %>% gather("tau", "missclass", -c("Sim", "Method", "n")) 

summary_peaks <- 
  peaks_long %>% group_by(Method, tau, n) %>% 
  summarise(
    mean_missclass = mean(missclass), 
    sd_missclass = sd(missclass)/sqrt(nSim)
  ) %>%
  ungroup() 


summary_peaks %>% 
  ggplot( aes(x = factor(n), y = mean_missclass, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_missclass - 2*sd_missclass, ymax = mean_missclass + 2*sd_missclass), 
                 position = position_dodge(width = 0.5))+
  facet_grid(.~factor(tau))+
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  ylim(c(0.09, .19)) +
  labs(x = "", y="Miss-classified rate",  col = "Method") 

ggsave("../Manuscript/Figures/peaks_missclass.png", width = 10, height = 3)


n <- 500
i <- 75
method <- "detrend_SIC"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             method, simDesign, n, i))
df$signal <- as.numeric(df$peaks > 0.5)
y_adj <- df$y - trend[,3]
signal_hat <- as.numeric(y_adj > thresh[3])
trend <- as.data.frame(trend)
trend$x <- df$x
ggplot(df, aes(x=x, y=y)) + geom_line() +
  geom_point(data=subset(df, signal==1), col="red", size = 2) +
  geom_point(data=df[signal_hat==1,], col="blue") + 
  geom_line(data=trend, aes(y=V3))
ggsave("../Manuscript/Figures/peaks_eg_class.png", width = 10, height = 3)

get_metric(trend[,3], df$y, df$signal, 
           thresh[3], "missclass")

get_metric(trend[,3], df$y, df$signal, 
           thresh[3], "precision")

sum(signal_hat==1&df$signal==1)/sum(signal_hat==1)

get_metric(trend[,3], df$y, df$signal, 
           thresh[4], "f1")
