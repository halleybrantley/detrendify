library(fields)
library(tidyverse)
library(devtools)
library(Hmisc)
library(caret)
library(zoo)
library(aricode)
library(mcclust)
load_all("detrendr")
rm(list=ls())
load("../SPod/trends_short.RData")
source("application_functions.R")
colPal <- c('#1b7837', '#c2a5cf')
            '#762a83','#9970ab',,
                '#a6dba0','#5aae61',)

spodPeaks <- select(spodPIDs, -time) - select(detrendr_trends,
                                              contains(paste(0.15)))

plot(spodPIDs$g, type="l")
lines(detrendr_trends$g_0.15, col="red")
lines(qsreg_trends$g_0.15, col="red")

plot(spodPeaks$h, type="l")
thresholds <- apply(spodPeaks, 2, 
                      function(x) median(x, na.rm=T) + 
                      2*median(abs(x-median(x, na.rm=T)), na.rm=T))
abline(h=thresholds[3], col="red")
plot(g~f, spodPeaks)

cor(spodPeaks$g, spodPeaks$f, method = "spearman")

################################################################################

methods <- c("detrendr", "qsreg")
metric_df <- tibble(metric=NA, method=NA, tau=NA, crit=NA, metric_type=NA)
i <- 1
metrics <- c("confusion", "NMI", "VI")
  
for (method in methods){
  trends <- get(paste(method, "trends", sep = "_"))
  for (j in 1:length(tau)){
    for (crit in c(2, 3, 4)){
      signal <- get_spod_signal(tau[j], trends, spodPIDs, crit)
      for (metric in metrics){
        if (metric == "confusion"){
          metric_df$metric[i] <- I(list(get_confusion(signal)))
        } else if (metric == "NMI"){
          metric_df$metric[i] <- I(list(get_NMI(signal)))
        } else if (metric == "VI") {
          metric_df$metric[i] <- I(list(get_VI(signal)))
        }
        metric_df$method[i] <- method
        metric_df$tau[i] <- tau[j]
        metric_df$crit[i] <- crit
        metric_df$metric_type[i] <- metric
        i <- i+1
        metric_df[i, ] <- NA
      }
    }
  }
}
metric_df <- metric_df[-i,]
  
NMI <- metric_df %>% 
  filter(metric_type %in% c("NMI")) %>% 
  select(-metric_type) %>% unnest() %>% 
  gather("nodes", "value", -c(method, tau, crit)) 

VI <- metric_df %>% 
  filter(metric_type %in% c("VI")) %>% 
  select(-metric_type) %>% unnest() %>% 
  gather("nodes", "value", -c(method, tau, crit)) 


ggplot(NMI, aes(x=nodes, y = value, col = method)) + 
  geom_point(position = position_dodge(width = .5), size = 2) +
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  facet_grid(crit~factor(tau)) +
  labs(y="Normalized Mutual Information", x = "Sensor Nodes", col = "")
ggsave("../Manuscript/Figures/NMI_app_short.png", width = 7, height = 3.5)



ggplot(VI, aes(x=nodes, y = value, col = method)) + 
  geom_point(position = position_dodge(width = .5), size = 2) +
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  facet_grid(crit~factor(tau)) +
  labs(y="Variation of Information", x = "Sensor Nodes", col = "")
ggsave("../Manuscript/Figures/VI_app_short.png", width = 7, height = 3.5)

################################################################################

confusion <- metric_df %>% 
  filter(metric_type %in% c("confusion")) %>% 
  select(-metric_type) %>% unnest() %>% 
  arrange(tau, crit)
  

latex(confusion %>% filter(crit == 3) %>% select(-crit),
      file = "../Manuscript/short_confusion_detrend_MAD3.tex",
      rowlabel = "",
      rowname = "",
      colheads = c("Method", "Quantile", 
                   "0,0,0", "1,0,0", "0,1,0", "1,1,0", 
                   "1,0,0", "1,1,0", "1,0,1", "1,1,1"),
      caption = "Confusion matrices for 3 SPod nodes after baseline 
      removal (n=5000). Node order is f, g, h. The threshold for the signal was 
      set as the median + 3*MAD.")

latex(confusion %>% filter(crit == 4) %>% select(-crit),
      file = "../Manuscript/short_confusion_detrend_MAD4.tex",
      rowlabel = "",
      rowname = "",
      colheads = c("Method", "Quantile", 
                   "0,0,0", "1,0,0", "0,1,0", "1,1,0", 
                   "1,0,0", "1,1,0", "1,0,1", "1,1,1"),
      caption = "Confusion matrices for 3 SPod nodes after baseline 
      removal (n=5000). Node order is f, g, h. The threshold for the signal was 
      set as the median + 4*MAD.")

latex(confusion %>% filter(crit == 2) %>% select(-crit),
      file = "../Manuscript/short_confusion_detrend_MAD5.tex",
      rowlabel = "",
      rowname = "",
      colheads = c("Method", "Quantile", 
                   "0,0,0", "1,0,0", "0,1,0", "1,1,0", 
                   "1,0,0", "1,1,0", "1,0,1", "1,1,1"),
      caption = "Confusion matrices for 3 SPod nodes after baseline 
      removal (n=5000). Node order is f, g, h. The threshold for the signal was 
      set as the median + 2*MAD.")

################################################################################
# Rug plot

plot(spodPIDs$f, type="l")
lines(qsreg_trends$f_0.15, col="red")
lines(detrendr_trends$f_0.15, col="red")
tmp1 <- spodPIDs$f - qsreg_trends$f_0.15
tmp2 <- spodPIDs$f - detrendr_trends$f_0.15

plot(tmp2, type="l")
lines(tmp1, col="red")

spod_signal <- get_spod_signal(0.15, qsreg_trends, spodPIDs, crit = 2)
spod_signal$time <- spodPIDs$time
spod_signal <- spod_signal %>%  gather(node, PID, -time) %>% filter(!is.na(node))
spodPeaks$time <- spodPIDs$time
spodLong <- spodPeaks %>% gather("node","PID", -time) %>% filter(!is.na(node))

spodLong$node <- factor(spodLong$node)
spodLong_signal <- spodLong[which(spod_signal$PID==1), ]

ggplot(spodLong, aes(x=time, y=PID)) +
  geom_rug(data = spodLong_signal, sides = "b") +
  geom_line(col="grey") +
  theme_bw() +
  facet_grid(node~., scales="free_y")
ggsave("../Manuscript/Figures/corrected_rugplot.png", width = 7, height = 3)
