################################################################################
# Calculate mis-classification rates in peaks simulation
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
rm(list=ls())
load_all("detrendr")
colPal <- rev(c('#762a83','#9970ab','#c2a5cf', '#d9f0d3',
                '#a6dba0','#5aae61','#1b7837'))

tau <- c(0.01, 0.05, 0.1)
nSim <- 100
methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "detrend_Xing", 
             "rqss", "npqw", "qsreg") 
metrics <- data.frame(method = NA, n = NA, sim = NA, 
                      metric = NA, metric_type = NA, threshold = NA)

thresh <- c(.9, 1, 1.1, 1.2)
metric_types <- "CAA" #c("missclass", "precision", "recall", "CAA", "F1")
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
            metrics$threshold[k] <- threshold
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
    save(metrics, file = "../SimResults/peak_metrics.RData")
  }
}
################################################################################
# Need to update

load("../SimResults/peak_metrics.RData")

noSignal <- metrics %>% 
  filter(metric_type == "recall", is.na(metric)) %>% 
  select(n, sim) %>% 
  distinct()

metrics$metric[which(is.na(metrics$metric))] <- 0

summary_peaks <- 
  metrics %>% 
  filter(!is.na(method)) %>%
  group_by(method, tau, n, threshold, metric_type) %>% 
  summarise(
    n_metrics = n(),
    mean_metric = mean(metric, na.rm=T), 
    sd_metric = sd(metric, na.rm=T)/sqrt(n_metrics)
  ) %>%
  ungroup() %>%
  mutate(method = factor(method, levels = methods))



for (metric_choice in metric_types){
  summary_peaks %>% 
    filter(metric_type == metric_choice) %>%
    ggplot( aes(x = factor(n), y = mean_metric, col = method)) + 
    geom_point(position = position_dodge(width = 0.5)) +
    geom_linerange(aes(ymin = mean_metric - 2*sd_metric, ymax = mean_metric + 
                         2*sd_metric), 
                   position = position_dodge(width = 0.5))+
    facet_grid(factor(threshold)~factor(tau))+
    theme_bw() +
    scale_color_manual(values=colPal, breaks = methods) +
    labs(x = "", y=metric_choice,  col = "Method") 
  
  ggsave(sprintf("../Manuscript/Figures/peaks_%s.png", metric_choice), 
         width = 7, height = 4.5)
}

n <- 1000
i <- 5
method <- "detrend_eBIC"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", "peaks", n, i))
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             method, "peaks", n, i))
df$signal <- as.numeric(df$peaks > 0.5)
y_adj <- df$y - trend[,2]
signal_hat <- as.numeric(y_adj > thresh[4])
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

get_metric(trend[,3], df$y, df$signal, 
           thresh[3], "missclass")

get_metric(trend[,3], df$y, df$signal, 
           thresh[3], "precision")

sum(signal_hat==1&df$signal==1)/sum(signal_hat==1)

get_metric(trend[,3], df$y, df$signal, 
           thresh[4], "f1")
