################################################################################
# Calculate mis-classification rates in peaks simulation
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
load_all("detrendr")
rm(list=ls())
spod <- read.csv("../SPod/fhrdata_2017-11-30.csv", 
                 header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

window_size <- 5000
overlap <- 500
max_iter <- 3
tau <- c(0.1, 0.2, 0.3)
k <- 3
spod_trends <- data.frame(time = spod$time)

for (node in c("f", "g", "h")){
  pidCol <- paste(node, "SPOD.PID..V.", sep=".")
  spodNode <- spod[, c("time", pidCol)]
  names(spodNode)[2] <- c("pid")
  spodNode$pid <- as.numeric(scale(spodNode$pid, center = TRUE))
  result <- get_windows_BIC(spodNode$pid, tau, k, window_size, overlap,
                          lambdaSeq = window_size^seq(1.1, 1.5, length.out=10), 
                          df_tol = 1e-9, 
                          gamma = 1,
                          plot_lambda = FALSE, 
                          solver = NULL, 
                          criteria = "eBIC")
  save(result, file=sprintf("../SPod/node_%s_trend.RData", node))
  spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
  names(spod_trends)[(ncol(spod_trends)-2):ncol(spod_trends)] <- 
    paste(node, tau, sep = "_")
}
save(spod_trends, file = "../SPod/spod_trends.RData")
load("../SPod/spod_trends.RData")
nodes <- c("f", "g", "h")
spodPIDs <- as.data.frame(scale(spod[, paste(nodes, "SPOD.PID..V.", sep=".")], 
                                center=TRUE))
names(spodPIDs) <- nodes
spodPeaks <- spodPIDs - select(spod_trends, contains("0.1"))

spodPeaks$time <- spod_trends$time  
#tmp <- movavg(na.locf(spodPeaks$f), 60, type="e")
#plot(spodPeaks$f[34500:36000], type="l")
#lines(tmp[34530:36030], col="blue")
#lines(spodPeaks$g[34500:36000], col="darkgreen")

spodPIDs$time <- spod_trends$time

spodRaw <- spodPIDs %>% gather("node", "PID", -time)
ggplot(spodRaw, aes(x=time, y=PID, col=node)) + 
  geom_line(alpha=0.5) + 
  theme_bw() 
ggsave("../Manuscript/Figures/uncorrected_data.png", width = 7, height = 3)

spodLong <- spodPeaks %>% gather("node","PID", -time)
thresholds <- 4*apply(spodPeaks[1:10000,1:3], 2, sd) + 
  apply(spodPeaks[1:10000,1:3], 2, mean)
PID_thresh <- data.frame(node = names(thresholds), thresh = thresholds)

ggplot(spodLong, aes(x=time, y=PID, col=node)) + 
  geom_line(alpha=0.5) + 
  geom_hline(data=PID_thresh, aes(yintercept=thresh, col=node)) +
  theme_bw() 
ggsave("../Manuscript/Figures/corrected_data.png", width = 7, height = 3)

ggplot(spodLong, aes(x=time, y=PID, col=node)) + 
  geom_line(alpha=0.5) + 
  theme_bw() +
  xlim(c(as.POSIXct("2017-11-30 9:09:57"), 
         as.POSIXct("2017-11-30 12:29:57")))
ggsave("../Manuscript/Figures/corrected_zoom_data.png", width = 7, height = 3)

head(spodPeaks)

get_spod_signal <- function(tau, spod_trends, spodPIDs){
  spodPeaks <- select(spodPIDs, -time) - select(spod_trends, contains(paste(tau)))
  thresholds <- 3*apply(spodPeaks[1:15000,1:3], 2, sd) + 
    apply(spodPeaks[1:15000,1:3], 2, mean)
  spodSignal <- spodPeaks
  for (i in 1:length(thresholds)){
    spodSignal[,i] <- as.numeric(spodPeaks[,i]>thresholds[i])
  }
  return(spodSignal)
}

get_misclass <- function(spodSignal){
  df <- data.frame(Comp = c("fg", "fh", "gh"), 
                   Missclass = 
                     c(mean(abs(spodSignal$f - spodSignal$g), na.rm=T),
                      mean(abs(spodSignal$f - spodSignal$h), na.rm=T),
                      mean(abs(spodSignal$g - spodSignal$h), na.rm=T)))
 return(df)
}

get_signal_ct <- function(spodSignal){
  df <- data.frame(Comp = c("fh", "gh", "fg", "fgh"), 
                   signal = 
                     c(sum(spodSignal$f == 1 & spodSignal$h == 1),
                       sum(spodSignal$g == 1 & spodSignal$h == 1), 
                       sum(spodSignal$f == 1 & spodSignal$g == 1), 
                       sum(spodSignal$f == 1 & spodSignal$g == 1 & spodSignal$h == 1)))
  return(df)
}


miss_df <- 
  cbind(get_misclass(get_spod_signal(0.1, spod_trends, spodPIDs)),
  get_misclass(get_spod_signal(0.2, spod_trends, spodPIDs))[,2],
  get_misclass(get_spod_signal(0.3, spod_trends, spodPIDs))[,2])

signal_comb <- 
  cbind(get_signal_ct(get_spod_signal(0.1, spod_trends, spodPIDs)),
        get_signal_ct(get_spod_signal(0.2, spod_trends, spodPIDs))[,2],
        get_signal_ct(get_spod_signal(0.3, spod_trends, spodPIDs))[,2])
names(signal_comb) <- c("Node", "tau 0.1", "tau 0.2", "tau 0.3")


signal_ct <- 
  tibble(Node = c("f", "g", "h"),
    `tau 0.1` = colSums(get_spod_signal(0.1, spod_trends, spodPIDs), na.rm=T),
    `tau 0.2` = colSums(get_spod_signal(0.2, spod_trends, spodPIDs), na.rm=T),
    `tau 0.3` = colSums(get_spod_signal(0.3, spod_trends, spodPIDs), na.rm=T))

latex(bind_rows(signal_ct, signal_comb), rowname="", title = '', 
      file = sprintf("../Manuscript/application_signal.tex"), 
      caption = "Seconds of signal by node combination and quantile level.")

names(miss_df) <- c("Comparison", "tau 0.1", "tau 0.2", "tau 0.3")
miss_df[,2:4] <- round(miss_df[,2:4],3)
latex(miss_df, rowname="", title = '', 
        file = sprintf("../Manuscript/application_missclass.tex"), 
      caption = "Fraction of seconds with different signal (0/1) classifications.")
