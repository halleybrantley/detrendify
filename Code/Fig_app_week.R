

library(tidyverse)
library(devtools)
library(Hmisc)
rm(list=ls())
source("application_functions.R")
colPal <- c('#1b7837', '#762a83')
VI <- {}
tau <- c(0.01, 0.05, 0.1)
nodes <- c("c", "d", "e")
cor.df <- data.frame(d = NA, tau = NA, 
                     cor_detrend = NA, cor_qsreg = NA)
class.df <- data.frame(d = NA, tau = NA, qthresh = NA,
                       VI_detrend = NA, VI_qsreg = NA, 
                       conf_detrend = NA, conf_qsreg = NA)
k <- 1
p <- 1
  load(sprintf("../SPod/Results/trends_e_2017-04-13.RData"))
  load(sprintf("../SPod/Results/qsreg_trends_2017-04-13.RData"))
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
}
class.df <- class.df[-p,] 
cor.df <- cor.df[-k,] 

VI <- class.df %>% 
  select(-c(d,starts_with("conf"))) %>% 
  gather(method, value, -c(tau, qthresh)) %>%
  unnest() %>% 
  gather("nodes", "value", -c(method, tau, qthresh)) 

ggplot(VI, aes(x=nodes, y = value, col = method)) + 
  geom_point(position = position_dodge(width = .5), size = 2) +
  theme_bw() +
  scale_color_manual(values=colPal, breaks = c("VI_detrend", "VI_qsreg")) +
  scale_x_discrete(labels = c("ab", "ac", "bc"))+
  facet_grid(qthresh~factor(tau), scales = "free") +
  labs(y="Variation of Information", x = "Sensor Nodes", col = "")

ggsave("../Manuscript/Figures/VI_by_day.png", width = 6, height = 5)    

# plot(c~e, peaks_qsreg)
# abline(v=thresh_1, col="red")
# abline(h=get_thresh2(thresh_1, coef_qsreg), col="red")
# 
plot(c~e, peaks_detrend)
# abline(v=thresh_1, col="red")
# abline(h=get_thresh2(thresh_1, coef_detrend), col="red")

VI$perc_signal <- VI$mean_pos/7200

ggplot(VI, aes(x=perc_signal, y=detrend)) +
  geom_point(aes(col="detrend"))+
  geom_point(aes(y=qsreg, col="qsreg")) +
  labs(y = "VI", x = "Proportion Signal")
