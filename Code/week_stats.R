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
spod_files <- dir("../SPod/SPod_week", ".RData", full.names=T)
colPal <- c('#1b7837', '#762a83')
nodes <- c("c", "e")
tau <- c(0.1, 0.15)
metric_all <- {}
for (i in 1:length(spod_files)){
  print(i)
load(spod_files[i])
# spodPeaks <- select(spodPIDs, -time) - select(spod_trends,
#                                               contains(paste(0.15)))
spodFig <- data.frame(time = spodPIDs$time,
                      detrendr = spod_trends$e_0.1,
                      qsreg = qsreg_trends$e_0.1) %>%
  gather("type","value", -time)

ggplot(spodFig, aes(x=time, y=value)) +
  #geom_line(data=spodPIDs, aes(y=c), col="darkgrey")+
  geom_line(aes(col=type, group=type)) +
  theme_bw() +
  scale_color_manual(breaks = c("detrendr", "qsreg"),
                     values = c(colPal)) +
  labs(col="", x="", y="PID") +
  xlim(c(spodPIDs$time[78000], spodPIDs$time[86400])) +
  ylim(c(1.2, 1.6))
plot(result$BIC[,1])

################################################################################
detrendr_trends <- spod_trends 
methods <- c("detrendr", "qsreg")
metric_df <- tibble(metric=NA, method=NA, tau=NA, crit=NA, metric_type=NA, h=NA)
i <- 1
metrics <- c("confusion",  "VI")
for (h in 0:23){
  start_ind <- h*3600+1
  end_ind <- min((h+1)*3600, nrow(spodPIDs))
  for (method in methods){
    trends <- get(paste(method, "trends", sep = "_"))[start_ind:end_ind, ]
    for (j in 1:length(tau)){
      for (crit in c(3, 4, 5)){
        signal <- get_spod_signal(tau[j], trends, 
                                  spodPIDs[start_ind:end_ind, ], 
                                  crit)
        for (metric in metrics){
          if (metric == "confusion"){
            metric_df$metric[i] <- I(list(get_confusion(signal, nodes)))
          } else if (metric == "VI") {
            metric_df$metric[i] <- as.numeric(get_VI(signal, nodes))
          }
          metric_df$method[i] <- method
          metric_df$tau[i] <- tau[j]
          metric_df$crit[i] <- crit
          metric_df$metric_type[i] <- metric
          metric_df$h = h
          i <- i+1
          metric_df[i, ] <- NA
        }
      }
    }
  }
}
metric_df <- metric_df[-i,]
metric_all <- rbind(metric_all, metric_df)
}

spodPeaks <- select(spodPIDs[start_ind:end_ind,], -time) - 
  select(spod_trends[start_ind:end_ind,], ends_with(paste(tau[1])))

par(mfrow=c(2,1))
plot(spodPIDs[start_ind:end_ind,"c"], type="l")
lines(spod_trends[start_ind:end_ind,"c_0.1"], col="red")
lines(qsreg_trends[start_ind:end_ind,"c_0.1"], col="blue")
plot(spodPIDs[start_ind:end_ind,"e"], type="l")
lines(spod_trends[start_ind:end_ind,"e_0.1"], col="red")
lines(qsreg_trends[start_ind:end_ind,"e_0.1"], col="blue")

VI <- metric_df %>% 
  filter(metric_type %in% c("VI")) %>% 
  select(-metric_type) %>% unnest() 

ggplot(VI, aes(x=method, y = metric, col = method)) + 
  geom_boxplot(position = position_dodge(width = .5)) +
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  facet_grid(crit~factor(tau), scales = "free") +
  labs(y="Variation of Information", x = "Sensor Nodes", col = "")
ggsave("../Manuscript/Figures/VI_app_short.png", width = 7, height = 3.5)
