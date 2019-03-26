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

tau <- c(0.1, 0.15)
metric_all <- {}
for (d in 2:8){
  print(d)
  load(sprintf("../SPod/SPod_week/trends_2017-03-0%d.RData",d))
  load(sprintf("../SPod/SPod_week/qsreg_trends_2017-03-0%d.RData",d))
# spodPeaks <- select(spodPIDs, -time) - select(qsreg_trends,
#                                               contains(paste(0.15)))
spodFig <- data.frame(time = spodPIDs$time,
                      detrendr = spod_trends$e_0.1,
                      qsreg = qsreg_trends$e_0.1) %>%
  gather("type","value", -time)
# 
# spodFig2 <- data.frame(time = spodPIDs$time,
#                       detrendr_e = spodPIDs$e - spod_trends$e_0.1,
#                       qsreg_e = spodPIDs$e - qsreg_trends$e_0.1, 
#                       detrendr_c = spodPIDs$c - spod_trends$c_0.1,
#                       qsreg_c = spodPIDs$c - qsreg_trends$c_0.1) 
#   gather("type","value", -time)
# 
# 
ggplot(spodFig, aes(x=time, y=value)) +
  geom_line(data=spodPIDs, aes(y=c), col="darkgrey")+
  geom_line(aes(col=type, group=type)) +
  theme_bw() +
  scale_color_manual(breaks = c("detrendr", "qsreg"),
                     values = c(colPal)) +
  labs(col="", x="", y="PID") 
  xlim(c(spodPIDs$time[1], spodPIDs$time[10000])) +
  ylim(c(1, 2))
# 
# plot(spodPeaks$c[80000:86400], type="l")
# lines(spodPeaks$e[80000:86400], col="red")
################################################################################
  detrendr_trends <- spod_trends 
  methods <- c("detrendr", "qsreg")
  metric_df <- tibble(metric=NA, method=NA, tau=NA, crit=NA, metric_type=NA, 
                      h=NA, d=NA)
  i <- 1
  metrics <- c("confusion",  "VI", "pos")
for (h in 0:2){
  start_ind <- h*28800+1
  end_ind <- min((h+1)*28800, nrow(spodPIDs))
  for (method in methods){
    trends <- get(paste(method, "trends", sep = "_"))[start_ind:end_ind, ]
    for (j in 1:length(tau)){
      for (crit in c(3, 4, 5)){
        signal <- get_spod_signal(tau[j], trends, 
                                  spodPIDs[start_ind:end_ind, ], 
                                  crit)
        
        for (metric in metrics){
          if (any(signal > 0)){
            if (metric == "confusion"){
              metric_df$metric[i] <- I(list(get_confusion(signal, nodes)))
            } else if (metric == "VI") {
              metric_df$metric[i] <- as.numeric(get_VI(signal, nodes))
            } else if (metric == "pos"){
              metric_df$metric[i] <- sum(rowSums(signal)>0)
            }
          }
          metric_df$method[i] <- method
          metric_df$tau[i] <- tau[j]
          metric_df$crit[i] <- crit
          metric_df$metric_type[i] <- metric
          metric_df$h[i] <- h
          metric_df$d[i] <- d
          i <- i+1
          metric_df[i, ] <- NA
        }
      }
    }
  }
}
metric_df <- metric_df[-i,]
metric_all <- bind_rows(metric_all, metric_df)
}

spodPeaks <- select(spodPIDs[start_ind:end_ind,], -time) -
  select(spod_trends[start_ind:end_ind,], ends_with(paste(tau[1])))

spodPeaks2 <- select(spodPIDs[start_ind:end_ind,], -time) -
  select(qsreg_trends[start_ind:end_ind,], ends_with(paste(tau[1])))
thresholds <- apply(spodPeaks, 2, 
                    function(x) median(x, na.rm=T) + 
                      crit*median(abs(x-median(x, na.rm=T)), na.rm=T))

par(mfrow=c(2,1))
plot(spodPeaks$c[1:15000], type="l")
abline(h=thresholds[1])
plot(spodPeaks$e[1:15000], type="l")
abline(h=thresholds[2])

thresholds2 <- apply(spodPeaks2, 2, 
                    function(x) median(x, na.rm=T) + 
                      crit*median(abs(x-median(x, na.rm=T)), na.rm=T))
par(mfrow=c(2,1))
plot(spodPeaks2$c[1:15000], type="l")
abline(h=thresholds2[1])
plot(spodPeaks2$e[1:15000], type="l")
abline(h=thresholds2[2])


par(mfrow=c(2,1))
plot(spodPIDs[start_ind:end_ind,"c"], type="l")
lines(spod_trends[start_ind:end_ind,"c_0.1"], col="red")
lines(qsreg_trends[start_ind:end_ind,"c_0.1"], col="blue")
plot(spodPIDs[start_ind:end_ind,"e"], type="l")
lines(spod_trends[start_ind:end_ind,"e_0.1"], col="red")
lines(qsreg_trends[start_ind:end_ind,"e_0.1"], col="blue")

VI <- metric_all %>% 
  filter(metric_type %in% c("VI")) %>% unnest()

VI %>% filter(d==3)%>%
ggplot(aes(x=h, y = metric, col = method)) + 
  geom_point() +
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  facet_grid(crit~factor(tau), scales = "free") +
  labs(y="Variation of Information", x = "Sensor Nodes", col = "")
ggsave("../Manuscript/Figures/VI_app_short.png", width = 7, height = 3.5)
