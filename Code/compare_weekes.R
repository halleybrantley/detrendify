

library(tidyverse)
library(devtools)
library(Hmisc)
rm(list=ls())
source("application_functions.R")
VI <- {}
for (d in 4:8){
load(sprintf("../SPod/Results/trends_e_2017-03-0%d.RData", d))
load(sprintf("../SPod/Results/qsreg_trends_2017-03-0%d.RData", d))

peaks_qsreg <- select(spodPIDs, -time) - 
  select(qsreg_trends, ends_with(paste(0.05)))

peaks_detrend <- select(spodPIDs, -time) - 
  select(spod_trends, ends_with(paste(0.05)))

# plot(c~e, peaks_qsreg)
# plot(c~e, peaks_detrend)

nodes <- c("c", "e")
cor(peaks_detrend[,nodes], method="spearman")
cor(peaks_qsreg[,nodes], method="spearman")

signal_detrend <- peaks_detrend
signal_qsreg <- peaks_qsreg
for (node in nodes){
  signal_detrend[,node] <- as.numeric(peaks_detrend[,node]>
                                        getCutoff(peaks_detrend, node))
  signal_qsreg[,node] <- as.numeric(peaks_qsreg[,node]>
                                      getCutoff(peaks_qsreg, node))
}

get_VI(signal_detrend, nodes)
get_VI(signal_qsreg, nodes)

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
  geom_point(aes(y=qsreg, col="qsreg")) 

plot(spodPIDs$c[1:ind_end], type="l", ylim=c(0,5))
lines(qsreg_trends$c_0.05[ind_start:ind_end], col="blue")
lines(spod_trends$c_0.05[ind_start:ind_end], col="red")

lines(spodPIDs$e[1:ind_end], type="l", col='red')
lines(qsreg_trends$e_0.05[ind_start:ind_end], col="blue")
lines(spod_trends$e_0.05[ind_start:ind_end], col="red")

thresh_e <- getCutoff(peaks_detrend, "e")
plot(peaks_detrend$e[ind_start:ind_end], type="l")
abline(h=thresh_e, col="red")

thresh_c <- getCutoff(peaks_detrend, "c")
lines(peaks_detrend$c[ind_start:ind_end], type="l", col="blue")
abline(h=thresh_c, col="red")
