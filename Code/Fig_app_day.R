################################################################################
# Figures for full day of data
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
rm(list=ls())
source("application_functions.R")
colPal <- c('#1b7837', '#762a83')
tau <- c(0.01, 0.05, 0.1)
nodes <- c("c", "d", "e")

load(sprintf("../SPod/Results/trends_e_2017-04-13.RData"))
load(sprintf("../SPod/Results/qsreg_trends_2017-04-13.RData"))

spodRaw <- spodPIDs %>%
  gather("node", "PID", -time)

ggplot(spodRaw, aes(x=time, y=PID, col=node)) +
  geom_line(alpha=0.5) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", labels = c("a", "b", "c"))+
  labs(col="SPod")

ggsave("../Manuscript/Figures/uncorrected_data.png", width = 7, height = 2.5)
################################################################################
peaks_qsreg <- select(spodPIDs, -time) - 
  select(qsreg_trends, ends_with(paste(tau[2])))

peaks_detrend <- select(spodPIDs, -time) - 
  select(spod_trends, ends_with(paste(tau[2])))

ggplot(spodPIDs, aes(x=c, y=e)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  labs(x = "SPod a", y = "SPod c", title = "Raw Data")
ggsave("../Manuscript/Figures/scatter_day_raw.png", width = 3, height = 3)

ggplot(peaks_detrend, aes(x=c, y=e)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  labs(x = "SPod a", y = "SPod c", title = "detrendr")
ggsave("../Manuscript/Figures/scatter_day_detrendr.png", width = 3, height = 3)

ggplot(peaks_qsreg, aes(x=c, y=e)) + 
  geom_point(alpha = 0.1) + 
  theme_bw() + 
  labs(x = "SPod a", y = "SPod c", title = "qsreg")
ggsave("../Manuscript/Figures/scatter_day_qsreg.png", width = 3, height = 3)
################################################################################  
VI <- {}
cor.df <- data.frame(tau = NA, 
                     cor_detrend = NA, cor_qsreg = NA)
class.df <- data.frame(tau = NA, qthresh = NA,
                       VI_detrend = NA, VI_qsreg = NA, 
                       conf_detrend = NA, conf_qsreg = NA)
k <- 1
p <- 1

for (j in 1:length(tau)){
  
  peaks_qsreg <- select(spodPIDs, -time) - 
    select(qsreg_trends, ends_with(paste(tau[j])))
  
  peaks_detrend <- select(spodPIDs, -time) - 
    select(spod_trends, ends_with(paste(tau[j])))

  cor.df$tau[k] <- tau[j]
  cor.df$cor_detrend[k] <- list(cor(peaks_detrend[,nodes], method="spearman",
                               use = "pairwise.complete.obs")[c(2,3,6)])
  cor.df$cor_qsreg[k] <- list(cor(peaks_qsreg[,nodes], method="spearman", 
                             use = "pairwise.complete.obs")[c(2,3,6)])
  k <- k+1
  cor.df[k,] <- NA
  
  for (q in c(0.9, 0.95, 0.99)){
    class.df$tau[p] <- tau[j]
    class.df$qthresh[p] <- q
    
    signal_detrend <- peaks_detrend
    signal_qsreg <- peaks_qsreg
    for (node in nodes){
      signal_detrend[,node] <- as.numeric(peaks_detrend[,node]>
                                            quantile(peaks_detrend[,node], 
                                                     q, na.rm=TRUE))
      signal_qsreg[,node] <- as.numeric(peaks_qsreg[,node]>
                                          quantile(peaks_qsreg[,node], 
                                                   q, na.rm=TRUE))
    }
    
    class.df$VI_detrend[p] <- list(get_VI(signal_detrend, nodes))
    class.df$VI_qsreg[p] <- list(get_VI(signal_qsreg, nodes))
    class.df$conf_detrend[p] <- list(get_confusion(signal_detrend, nodes))
    class.df$conf_qsreg[p] <- list(get_confusion(signal_qsreg, nodes))
    p <- p+1
    class.df[p,] <- NA
  }
}

class.df <- class.df[-p,] 
cor.df <- cor.df[-k,] 

VI <- class.df %>% 
  select(-c(starts_with("conf"))) %>% 
  gather(method, value, -c(tau, qthresh)) %>%
  unnest() %>% 
  gather("nodes", "value", -c(method, tau, qthresh)) 

ggplot(VI, aes(x=nodes, y = value, col = method, shape = factor(qthresh))) + 
  geom_point(position = position_dodge(width = .5), size = 2) +
  theme_bw() +
  scale_color_manual(values=colPal, breaks =c("VI_detrend", "VI_qsreg"), 
                      labels = c("detrendr", "qsreg")) +
  scale_x_discrete(labels = c("ab", "ac", "bc"))+
  facet_grid(~factor(tau)) +
  labs(y="Variation of Information", x = "SPods", col = "Method", 
       shape = "Threshold Quantile")
ggsave("../Manuscript/Figures/VI_app_long.png", width = 7, height = 3.5) 

################################################################################

