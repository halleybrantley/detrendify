library(fields)
library(tidyverse)
library(devtools)
library(Hmisc)
library(caret)
library(aricode)
library(mcclust)
library(lubridate)
load_all("detrendr")
rm(list=ls())
source("application_functions.R")
colPal <- c('#1b7837', '#762a83')


load(sprintf("../SPod/SPod_week/trends_e_2017-04-13.RData"))
load(sprintf("../SPod/SPod_week/qsreg_trends_2017-04-13.RData"))

# spodFig <- data.frame(time = spodPIDs$time,
#                      # detrendr = spod_trends$e_0.1,
#                       qsreg = qsreg_trends$e_0.01) %>%
#   gather("type","value", -time)
# 
# ggplot(spodFig, aes(x=time, y=value)) +
#   geom_line(data=spodPIDs, aes(y=e), col="darkgrey")+
#   geom_line(aes(col=type, group=type)) +
#   theme_bw() +
#   scale_color_manual(breaks = c("detrendr", "qsreg"),
#                      values = c(colPal)) +
#   labs(col="", x="", y="PID") +
# xlim(c(spodPIDs$time[1], spodPIDs$time[10000])) +
#   ylim(c(1, 2))

################################################################################
  detrendr_trends <- spod_trends 
  methods <- "qsreg" #c("detrendr", "qsreg")
  metric_df <- tibble(metric=NA, method=NA, tau=NA, crit=NA, metric_type=NA, 
                      h=NA)
  i <- 1
  metrics <- c("VI", "confusion")
  for (h in 0:1){
    start_ind <- h*43200+1
    end_ind <- min((h+1)*43200, nrow(spodPIDs))
    for (method in methods){
      trends <- get(paste(method, "trends", sep = "_"))[start_ind:end_ind, ]
      for (j in 1:length(tau)){
        for (crit in c(5, 6)){
          signal <- get_spod_signal(tau[j], trends, 
                                    spodPIDs[start_ind:end_ind, ], 
                                    crit)
          
          for (metric in metrics){
            if (any(signal > 0)){
              if (metric == "confusion"){
                metric_df$metric[i] <- I(list(get_confusion(signal, nodes)))
              } else if (metric == "VI") {
                metric_df$metric[i] <- I(list(get_VI(signal, nodes)))
              } else if (metric == "pos"){
                metric_df$metric[i] <- sum(rowSums(signal)>0)
              }
            }
            metric_df$method[i] <- method
            metric_df$tau[i] <- tau[j]
            metric_df$crit[i] <- crit
            metric_df$metric_type[i] <- metric
            metric_df$h[i] <- h
            i <- i+1
            metric_df[i, ] <- NA
          }
        }
      }
    }
  }

VI <- metric_df %>% 
  filter(metric_type %in% c("VI")) %>% 
  select(-metric_type) %>% unnest() %>% 
  gather("nodes", "value", -c(method, tau, crit, h)) 

ggplot(VI, aes(x=nodes, y = value, col = method)) + 
  geom_point(position = position_dodge(width = .5), size = 2) +
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  scale_x_discrete(labels = c("ab", "ac", "bc"))+
  facet_grid(crit~factor(tau), scales = "free") +
  labs(y="Variation of Information", x = "Sensor Nodes", col = "")
ggsave("../Manuscript/Figures/VI_app_day.png", width = 7, height = 3.5)
