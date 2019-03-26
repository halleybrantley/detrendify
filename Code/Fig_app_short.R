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
load("../SPod/trends_short2.RData")
source("application_functions.R")
colPal <- c('#1b7837', '#762a83')
nodes <- c("c", "d", "e")

spodPeaks <- select(spodPIDs, -time) - select(detrendr_trends,
                                              contains(paste(0.15)))
spodFig <- data.frame(time = spodPIDs$time, 
                      detrendr = detrendr_trends$c_0.15, 
                      qsreg = qsreg_trends$c_0.15) %>%
  gather("type","value", -time)

ggplot(spodFig, aes(x=time, y=value)) + 
  geom_line(data=spodPIDs, aes(y=c), col="darkgrey")+
  geom_line(aes(col=type, group=type)) + 
  theme_bw() +
  scale_color_manual(breaks = c("detrendr", "qsreg"), 
                     values = c(colPal)) +
  labs(col="", x="", y="PID")
ggsave("../Manuscript/Figures/short_trends.png", width = 7, height = 2.5)


################################################################################

methods <- c("detrendr", "qsreg")
metric_df <- tibble(metric=NA, method=NA, tau=NA, crit=NA, metric_type=NA)
i <- 1
metrics <- c("confusion", "NMI", "VI")
  
for (method in methods){
  trends <- get(paste(method, "trends", sep = "_"))
  for (j in 1:length(tau)){
    for (crit in c(3, 4, 5)){
      signal <- get_spod_signal(tau[j], trends, spodPIDs, crit)
      for (metric in metrics){
        if (metric == "confusion"){
          metric_df$metric[i] <- I(list(get_confusion(signal, nodes)))
        } else if (metric == "NMI"){
          metric_df$metric[i] <- I(list(get_NMI(signal, nodes)))
        } else if (metric == "VI") {
          metric_df$metric[i] <- I(list(get_VI(signal, nodes)))
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
  scale_x_discrete(labels = c("ab", "ac", "bc"))+
  facet_grid(crit~factor(tau)) +
  labs(y="Normalized Mutual Information", x = "Sensor Nodes", col = "")
ggsave("../Manuscript/Figures/NMI_app_short.png", width = 7, height = 3.5)



ggplot(VI, aes(x=nodes, y = value, col = method)) + 
  geom_point(position = position_dodge(width = .5), size = 2) +
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  scale_x_discrete(labels = c("ab", "ac", "bc"))+
  facet_grid(crit~factor(tau), scales = "free") +
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
      removal (n=6000). Node order is a, b, c. The threshold for the signal was 
      set as the median + 3*MAD.")

latex(confusion %>% filter(crit == 4) %>% select(-crit),
      file = "../Manuscript/short_confusion_detrend_MAD4.tex",
      rowlabel = "",
      rowname = "",
      colheads = c("Method", "Quantile", 
                   "0,0,0", "1,0,0", "0,1,0", "1,1,0", 
                   "1,0,0", "1,1,0", "1,0,1", "1,1,1"),
      caption = "Confusion matrices for 3 SPod nodes after baseline 
      removal (n=6000). Node order is a, b, c. The threshold for the signal was 
      set as the median + 4*MAD.")

latex(confusion %>% filter(crit == 5) %>% select(-crit),
      file = "../Manuscript/short_confusion_detrend_MAD5.tex",
      rowlabel = "",
      rowname = "",
      colheads = c("Method", "Quantile", 
                   "0,0,0", "1,0,0", "0,1,0", "1,1,0", 
                   "1,0,0", "1,1,0", "1,0,1", "1,1,1"),
      caption = "Confusion matrices for 3 SPod nodes after baseline 
      removal (n=6000). Node order is a, b, c. The threshold for the signal was 
      set as the median + 5*MAD.")

################################################################################
# Rug plot
crit <- 5
tau0 <- 0.15
trends <- detrendr_trends

# par(mfrow=c(2,1), mar = c(2,2,0,0))
# plot(spodPIDs[, nodes[1]], type="l")
# lines(qsreg_trends$c_0.1, col="red")
# lines(detrendr_trends$c_0.1, col="red")
# plot(spodPIDs[, nodes[3]], type="l")
# lines(qsreg_trends$e_0.1, col="red")
# lines(detrendr_trends$e_0.1, col="red")
# tmp1 <- spodPIDs[, nodes[1]] - qsreg_trends[, "c_0.15"]
# tmp2 <- spodPIDs[, nodes[3]] - qsreg_trends[, "e_0.15"]
# 
# tmp3 <- spodPIDs[, nodes[1]] - detrendr_trends[, "c_0.15"]
# tmp4 <- spodPIDs[, nodes[3]] - detrendr_trends[, "e_0.15"]
# plot(tmp3, type="l")
# lines(tmp4, col="red")
# 
# plot(tmp3, type="l")
# lines(tmp4, col="red")
# cor(tmp3, tmp4, method = "spearman")
# cor(tmp1, tmp2, method = "spearman")

spod_signal <- get_spod_signal(tau0, trends, spodPIDs, crit = crit)
spod_signal$time <- spodPIDs$time
spod_signal <- spod_signal %>%  
  gather(node, PID, -time) %>% 
  filter(!is.na(node))
spodPeaks <- select(spodPIDs, -time) - 
  select(trends, contains(paste(tau0)))
thresholds <- data.frame(node = nodes, 
                         thresh = apply(spodPeaks, 2, 
                    function(x) median(x, na.rm=T) + 
                      crit*median(abs(x-median(x, na.rm=T)), na.rm=T)))
thresholds$node <- c("a", "b", "c")
spodPeaks$time <- spodPIDs$time
spodLong <- spodPeaks %>% 
  gather("node","PID", -time) %>% 
  filter(!is.na(node))

spodLong$node <- factor(spodLong$node, labels = c("a", "b", "c"))
spodLong_signal <- spodLong[which(spod_signal$PID==1), ]

ggplot(spodLong, aes(x=time, y=PID)) +
  geom_rug(data = spodLong_signal, sides = "b") +
  geom_line(col="grey") +
  theme_bw() +
  geom_hline(data=thresholds, aes(yintercept = thresh), linetype="dashed")+
  facet_grid(node~., scales="free_y") 
  # xlim(c(as.POSIXct("2017-11-30 11:00:01 EST"), 
  #        as.POSIXct("2017-11-30 12:00:01 EST"))) +
  #ylim(c(0, .2))
ggsave("../Manuscript/Figures/corrected_rugplot.png", width = 7, height = 4)
tail(spodLong[which(spodLong$PID > 0.2),])
