################################################################################
# Calculate mis-classification rates in peaks simulation and create figures
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(detrendr)
rm(list=ls())

colPal <- c('#006d2c', '#2ca25f', '#66c2a4', '#1c9099', '#d8b365',
            "#c2a5cf", "#9970ab", "#762a83")

methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "detrend_Xing", 
             "windows",
             "rqss", "npqw", "qsreg") 

tau <- c(0.01, 0.05, 0.1)
nSim <- 100
thresh <- c(.9, 1, 1.1, 1.2)
metric_types <- "CAA" #c("missclass", "precision", "recall", "CAA", "F1")

# metrics <- data.frame(method = NA, n = NA, sim = NA, 
#                       metric = NA, metric_type = NA, threshold = NA)
# k <- 1
# 
# for (n in c(500,1000,2000,4000)){
#   for (i in 1:nSim){
#     load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", "peaks", n, i))
#     df$signal <- as.numeric(df$peaks > 0.5)
# 
#     for (method in methods){
#       load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData",
#                    method, "peaks", n, i))
#       for (threshold in thresh){
#         for (metric in metric_types){
#           for (j in 1:length(tau)){
#             metrics$sim[k] <- i
#             metrics$method[k] <- method
#             metrics$n[k] <- n
#             metrics$threshold[k] <- threshold
#             metrics$tau[k] <- tau[j]
#             metrics$metric_type[k] <- metric
#             metrics$metric[k] <- get_metric(trend[,j], df$y, df$signal,
#                                          threshold, metric)
#             k <- k+1
#             metrics[k, ] <- NA
#           }
#         }
#       }
#     }
#     save(metrics, file = "../SimResults/peak_metrics.RData")
#   }
# }
# metrics <- metrics[-k,]
# save(metrics, file = "../SimResults/peak_metrics.RData")
################################################################################

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
  
  ggsave(sprintf("../Manuscript/Figures/peaks_with_windows_%s.png", metric_choice), 
         width = 7, height = 4.5)
  
  summary_peaks %>% 
    filter(metric_type == metric_choice, method != "windows") %>%
    ggplot( aes(x = factor(n), y = mean_metric, col = method)) + 
    geom_point(position = position_dodge(width = 0.5)) +
    geom_linerange(aes(ymin = mean_metric - 2*sd_metric, ymax = mean_metric + 
                         2*sd_metric), 
                   position = position_dodge(width = 0.5))+
    facet_grid(factor(threshold)~factor(tau))+
    theme_bw() +
    scale_color_manual(values=colPal[-5], breaks = methods) +
    labs(x = "", y=metric_choice,  col = "Method") 
  
  ggsave(sprintf("../Manuscript/Figures/peaks_%s.png", metric_choice), 
         width = 7, height = 4.5)
}


summary_stats <- summary_peaks %>% 
  filter(metric_type == metric_choice, method != "windows")


summary_stats$value <- sprintf("%0.3f (%0.3f)", summary_stats$mean_metric, 
                               summary_stats$sd_metric)

wide_stats <- 
  summary_stats %>% 
  select(method, tau, threshold, n, value) %>%
  spread(tau, value) %>%
  arrange(threshold, n, method)
wide_stats$method <- as.character(str_replace(wide_stats$method, "_", " ")) 
unique(wide_stats$method)

latex(wide_stats %>% filter(threshold==0.9) %>% select(-n, -threshold) , 
      file = "Fig10_09.tex",       
      rowname = "",
      title = '', 
      n.rgroup = c(7,7,7,7),
      rgroup = c("n=500", "n=1000", "n=2000", "n=4000"))


latex(wide_stats %>% filter(threshold==1.0) %>% select(-n, -threshold) , 
      file = "Fig10_10.tex",       
      rowname = "",
      title = '', 
      n.rgroup = c(7,7,7,7),
      rgroup = c("n=500", "n=1000", "n=2000", "n=4000"))

latex(wide_stats %>% filter(threshold==1.1) %>% select(-n, -threshold) , 
      file = "Fig10_11.tex",       
      rowname = "",
      title = '', 
      n.rgroup = c(7,7,7,7),
      rgroup = c("n=500", "n=1000", "n=2000", "n=4000"))

latex(wide_stats %>% filter(threshold==1.2) %>% select(-n, -threshold) , 
      file = "Fig10_12.tex",       
      rowname = "",
      title = '', 
      n.rgroup = c(7,7,7,7),
      rgroup = c("n=500", "n=1000", "n=2000", "n=4000"))
