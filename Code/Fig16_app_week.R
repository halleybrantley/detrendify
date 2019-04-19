################################################################################
# Figure for full week of data
# Halley Brantley
################################################################################

library(tidyverse)
library(devtools)
library(Hmisc)
library(caret)
rm(list=ls())
source("application_functions.R")
VI <- {}
tau <- c(0.01, 0.05, 0.1)
nodes <- c("c", "e")
cor.df <- data.frame(d = NA, tau = NA, 
                     cor_detrend = NA, cor_qsreg = NA)
class.df <- data.frame(d = NA, tau = NA, qthresh = NA,
                       VI_detrend = NA, VI_qsreg = NA, 
                       conf_detrend = NA, conf_qsreg = NA)
k <- 1
p <- 1
for (d in 2:8){
  load(sprintf("../SPod/Results/trends_j_2018-06-%d.RData", d+13))
  load(sprintf("../SPod/Results/qsreg_trends_2018-06-%d.RData", d+13))
  for (j in 1:length(tau)){

    peaks_qsreg <- select(spodPIDs, -time) - 
      select(qsreg_trends, ends_with(paste(tau[j])))
    
    peaks_detrend <- select(spodPIDs, -time) - 
      select(spod_trends, ends_with(paste(tau[j])))
    
    cor.df$d[k] <- d
    cor.df$tau[k] <- tau[j]
    cor.df$cor_detrend[k] <- cor(peaks_detrend[,nodes], method="spearman")[1,2]
    cor.df$cor_qsreg[k] <- cor(peaks_qsreg[,nodes], method="spearman")[1,2]
    k <- k+1
    cor.df[k,] <- NA
    
    for (q in c(0.9, 0.95, 0.99)){
      class.df$d[p] <- d
      class.df$tau[p] <- tau[j]
      class.df$qthresh[p] <- q
      
      signal_detrend <- peaks_detrend
      signal_qsreg <- peaks_qsreg
      for (node in nodes){
        signal_detrend[,node] <- as.numeric(peaks_detrend[,node]>
                                              quantile(peaks_detrend[,node], 
                                                       q))
        signal_qsreg[,node] <- as.numeric(peaks_qsreg[,node]>
                                            quantile(peaks_qsreg[,node], 
                                                     q))
      }

      class.df$VI_detrend[p] <- get_VI(signal_detrend, nodes)
      class.df$VI_qsreg[p] <- get_VI(signal_qsreg, nodes)
      class.df$conf_detrend[p] <- list(get_confusion(signal_detrend, nodes))
      class.df$conf_qsreg[p] <- list(get_confusion(signal_qsreg, nodes))
      p <- p+1
      class.df[p,] <- NA
    }
  }
}
class.df <- class.df[-p,] 
cor.df <- cor.df[-k,] 


ggplot(class.df, aes(x=VI_qsreg, y=VI_detrend, 
                     col=factor(qthresh), shape = factor(tau)))+
  geom_point(size = 2) + 
  geom_abline(slope=1, intercept = 0, linetype=2) +
  theme_bw() +
  scale_color_brewer(palette="Set1") +
  ylim(c(0,1)) +
  xlim(c(0,1)) +
  labs(x = "qsreg", y = "detrendr", col = "Quantile Threshold", 
       shape = "Quantile Baseline") 
ggsave("../Manuscript/Figures/VI_by_day.png", width = 5, height = 3.3)    







