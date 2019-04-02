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
colPal <- c('#1b7837', '#762a83', "darkgrey")
nodes <- c("c", "d", "e")
missInd <- which(is.na(spodPIDs$d))
tau <- c(0.01, 0.05, 0.1)
detrendr_trends[missInd, paste("d", tau, sep="_")] <- NA
qsreg_trends[missInd, paste("d", tau, sep="_")] <- NA
spodPeaks <- select(spodPIDs, -time) - select(detrendr_trends,
                                              contains(paste(0.05)))
spodFig <- data.frame(time = c(spodPIDs$time, spodPIDs$time, spodPIDs$time), 
                      raw = c(spodPIDs$c, spodPIDs$d, spodPIDs$e),
                      detrendr = c(detrendr_trends$c_0.05, 
                                   detrendr_trends$d_0.05, 
                                   detrendr_trends$e_0.05), 
                      qsreg = c(qsreg_trends$c_0.05, 
                                qsreg_trends$d_0.05, 
                                qsreg_trends$e_0.05), 
                      node = c(rep("a", nrow(spodPIDs)), 
                               rep("b", nrow(spodPIDs)), 
                               rep("c", nrow(spodPIDs)))) %>%
  gather("type","value", -c(time, node)) %>%
  arrange(time, rev(type))

ggplot(spodFig, aes(x=time, y=value)) + 
  geom_line(aes(col=type, group=type), alpha=0.5) + 
  theme_bw() +
  facet_grid(node~., scales="free")+
  scale_color_manual(breaks = c("raw", "detrendr", "qsreg"), 
                     values = c(colPal)) +
  labs(col="", x="", y="PID")
ggsave("../Manuscript/Figures/short_trends.png", width = 7, height = 4)


################################################################################
methods <- c("detrendr", "qsreg")
metric_df <- tibble(metric=NA, method=NA, tau=NA, crit=NA, metric_type=NA)
i <- 1
metrics <- c("confusion", "VI")
  
for (method in methods){
  trends <- get(paste(method, "trends", sep = "_"))
  for (j in 1:length(tau)){
    for (crit in c(0.90, 0.95, 0.99)){
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

VI <- metric_df %>% 
  filter(metric_type %in% c("VI")) %>% 
  select(-metric_type) %>% unnest() %>% 
  gather("nodes", "value", -c(method, tau, crit)) 

ggplot(VI, aes(x=nodes, y = value, col = method, shape = factor(crit))) + 
  geom_point(position = position_dodge(width = .5), size = 2) +
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  scale_x_discrete(labels = c("ab", "ac", "bc"))+
  facet_grid(~factor(tau)) +
  labs(y="Variation of Information", x = "SPods", col = "Method", 
       shape = "Threshold Quantile")
ggsave("../Manuscript/Figures/VI_app_short.png", width = 7, height = 3.5)
################################################################################

# Rug plot
crit <- .95
tau0 <- 0.05
trends <- detrendr_trends

spod_signal <- get_spod_signal(tau0, trends, spodPIDs, crit = crit)
spod_signal$time <- spodPIDs$time
spod_signal <- spod_signal %>%  
  gather(node, PID, -time) %>% 
  filter(!is.na(node))
spodPeaks <- select(spodPIDs, -time) - 
  select(trends, contains(paste(tau0)))
thresholds <- data.frame(node = nodes, 
                         thresh = apply(spodPeaks, 2, quantile, crit, na.rm=T))
thresholds$node <- c("a", "b", "c")
thresholds_90 <- data.frame(node = nodes, 
                            thresh = apply(spodPeaks, 2, quantile, 0.90, na.rm=T))
thresholds_99 <- data.frame(node = nodes, 
                            thresh = apply(spodPeaks, 2, quantile, 0.99, na.rm=T))
thresholds_90$node <- c("a", "b", "c")
thresholds_99$node <- c("a", "b", "c")

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
  geom_hline(data=thresholds_90, aes(yintercept = thresh), linetype="dashed", 
             col="blue") +
  geom_hline(data=thresholds_99, aes(yintercept = thresh), linetype="dashed", 
             col="orange") +
  facet_grid(node~., scales="free_y") 
  # xlim(c(as.POSIXct("2017-11-30 11:00:01 EST"), 
  #        as.POSIXct("2017-11-30 12:00:01 EST"))) +
  #ylim(c(0, .2))
ggsave("../Manuscript/Figures/corrected_rugplot.png", width = 7, height = 4)
tail(spodLong[which(spodLong$PID > 0.2),])

################################################################################

confusion <- metric_df %>% 
  filter(metric_type %in% c("confusion")) %>% 
  select(-metric_type) %>% unnest() %>% 
  arrange(tau, crit)


latex(confusion %>% filter(crit == 5) %>% select(-crit),
      file = "../Manuscript/short_confusion_detrend_MAD5.tex",
      rowlabel = "",
      rowname = "",
      colheads = c("Method", "Quantile", 
                   "0,0,0", "1,0,0", "0,1,0", "1,1,0", 
                   "1,0,0", "1,1,0", "1,0,1", "1,1,1"),
      caption = "Confusion matrices for 3 SPod nodes after baseline 
      removal (n=7200). Node order is a, b, c. The threshold for the signal was 
      set as the median + 5*MAD.")

latex(confusion %>% filter(crit == 6) %>% select(-crit),
      file = "../Manuscript/short_confusion_detrend_MAD6.tex",
      rowlabel = "",
      rowname = "",
      colheads = c("Method", "Quantile", 
                   "0,0,0", "1,0,0", "0,1,0", "1,1,0", 
                   "1,0,0", "1,1,0", "1,0,1", "1,1,1"),
      caption = "Confusion matrices for 3 SPod nodes after baseline 
      removal (n=7200). Node order is a, b, c. The threshold for the signal was 
      set as the median + 6*MAD.")

################################################################################
