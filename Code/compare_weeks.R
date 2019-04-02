

library(tidyverse)
library(devtools)
library(Hmisc)
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
for (d in 2:9){
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
  labs(x = "qsreg", y = "detrendr", col = "Quantile Threshold", 
       shape = "Quantile Baseline") 
ggsave("../Manuscript/Figures/VI_by_day.png", width = 5, height = 4)    

# plot(c~e, peaks_qsreg)
# abline(v=thresh_1, col="red")
# abline(h=get_thresh2(thresh_1, coef_qsreg), col="red")
# 
# plot(c~e, peaks_detrend)
# abline(v=thresh_1, col="red")
# abline(h=get_thresh2(thresh_1, coef_detrend), col="red")








VI.df <- data.frame(detrend = NA, qsreg = NA, mean_pos = NA)
for (j in 1:12){
  ind_start <- (j-1)*7200 + 1
  ind_end <- min(nrow(spodPIDs), j*7200)
  VI.df$mean_pos[j] <- sum(apply(cbind(signal_detrend[ind_start:ind_end,], 
                            signal_qsreg[ind_start:ind_end,]), 1, 
                      function(x) any(x>0)))
  VI.df$detrend[j] <- as.numeric(get_VI(signal_detrend[ind_start:ind_end,], 
                                        nodes))
  VI.df$qsreg[j] <- as.numeric(get_VI(signal_qsreg[ind_start:ind_end,], 
                                      nodes))
  VI.df[j+1, ] <- NA
}
VI.df <- VI.df[-(j+1), ]
VI <- rbind(VI, VI.df)
}

VI$perc_signal <- VI$mean_pos/7200

ggplot(VI, aes(x=perc_signal, y=detrend)) +
  geom_point(aes(col="detrend"))+
  geom_point(aes(y=qsreg, col="qsreg")) +
  labs(y = "VI", x = "Proportion Signal")
ind_start <- 7200
ind_end <- 7200*4
plot(spodPIDs[ind_start:ind_end, nodes[1]], type="l")
lines(qsreg_trends[ind_start:ind_end, paste0(nodes[1], "_0.05")], col="blue")
lines(spod_trends[ind_start:ind_end, paste0(nodes[1], "_0.05")], col="red")

plot(spodPIDs[ind_start:ind_end, nodes[2]], type="l", col='red')
lines(qsreg_trends[ind_start:ind_end, paste0(nodes[2], "_0.05")], col="blue")
lines(spod_trends[ind_start:ind_end, paste0(nodes[2], "_0.05")], col="red")

par(mfrow = c(2,1))
thresh_1 <- getCutoff(peaks_detrend, nodes[1])
plot(peaks_detrend[ind_start:ind_end, nodes[1]], type="l")
abline(h=thresh_1, col="red")

thresh_2 <- getCutoff(peaks_detrend, nodes[2])
plot(peaks_detrend[ind_start:ind_end, nodes[2]], type="l", col="blue")
abline(h=thresh_2, col="blue")

plot(peaks_detrend[, nodes[1]], type="l")
lines(peaks_detrend[, nodes[2]], type="l", col="red")
