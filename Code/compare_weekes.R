

library(tidyverse)
library(devtools)
library(Hmisc)
rm(list=ls())
source("application_functions.R")
VI <- {}
for (d in 4:8){
load(sprintf("../SPod/Results/trends_j_2018-06-%d.RData", d+13))
load(sprintf("../SPod/Results/qsreg_trends_2018-06-%d.RData", d+13))


peaks_qsreg <- select(spodPIDs, -time) - 
  select(qsreg_trends, ends_with(paste(0.05)))

peaks_detrend <- select(spodPIDs, -time) - 
  select(spod_trends, ends_with(paste(0.05)))

# plot(j~d, peaks_qsreg)
# abline(v=thresh_1, col="red")
# abline(h=get_thresh2(thresh_1, coef_qsreg), col="red")
# 
# plot(j~d, peaks_detrend)
# abline(v=thresh_1, col="red")
# abline(h=get_thresh2(thresh_1, coef_detrend), col="red")

nodes <- c("d", "j")
cor(peaks_detrend[,nodes], method="spearman")
cor(peaks_qsreg[,nodes], method="spearman")

signal_detrend <- peaks_detrend
signal_qsreg <- peaks_qsreg
qthresh <- 0.95
for (node in nodes){
  signal_detrend[,node] <- as.numeric(peaks_detrend[,node]>
                                          quantile(peaks_detrend[,node], 
                                                   qthresh))
  signal_qsreg[,node] <- as.numeric(peaks_qsreg[,node]>
                                        quantile(peaks_qsreg[,node], 
                                                 qthresh))
}


get_VI(signal_detrend, nodes)
get_confusion(signal_detrend, nodes)
get_VI(signal_qsreg, nodes)
get_confusion(signal_qsreg, nodes)

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
