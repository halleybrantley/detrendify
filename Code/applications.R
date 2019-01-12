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
# spod_trends <- data.frame(time = spod$time)

# for (node in c("f", "g", "h")){
#   pidCol <- paste(node, "SPOD.PID..V.", sep=".")
#   spodNode <- spod[, c("time", pidCol)]
#   names(spodNode)[2] <- c("pid")
#   spodNode$pid <- as.numeric(scale(spodNode$pid, center = TRUE))
#   result <- get_windows_BIC(spodNode$pid, tau, k, window_size, overlap,
#                           lambdaSeq = window_size^seq(1.1, 1.5, length.out=10), 
#                           df_tol = 1e-9, 
#                           gamma = 1,
#                           plot_lambda = FALSE, 
#                           solver = NULL, 
#                           criteria = "eBIC")
#   save(result, file=sprintf("../SPod/node_%s_trend.RData", node))
#   spod_trends <- cbind(spod_trends, as.data.frame(result$trend))
#   names(spod_trends)[(ncol(spod_trends)-2):ncol(spod_trends)] <- 
#     paste(node, tau, sep = "_")
# }
# save(spod_trends, file = "../SPod/spod_trends.RData")
################################################################################

load("../SPod/spod_trends.RData")
nodes <- c("f", "g", "h")
spodPIDs <- as.data.frame(scale(spod[, paste(nodes, "SPOD.PID..V.", sep=".")], 
                                center=TRUE))
names(spodPIDs) <- nodes
spodPeaks <- spodPIDs - select(spod_trends, contains("0.1"))
spodPeaks$time <- spod_trends$time  
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

################################################################################

get_spod_signal <- function(tau, spod_trends, spodPIDs, thresholds){
  spodPeaks <- select(spodPIDs, -time) - select(spod_trends, contains(paste(tau)))
  spodSignal <- spodPeaks
  for (i in 1:length(thresholds)){
    spodSignal[,i] <- as.numeric(spodPeaks[,i]>thresholds[i])
  }
  return(spodSignal)
}

spod_signal <- get_spod_signal(0.1, spod_trends, spodPIDs, PID_thresh$thresh)
get_metric(trend[,j], df$y, df$signal,
           threshold, metric)
metric_types <- c("CAA", "F1", "precision", "recall", "missclass")
metrics <- data.frame(metric = NA, metric_type = NA, tau = NA, 
                      node_true = NA, node_est = NA)

k <- 1
for (node1 in nodes){
  for (node2 in nodes){
    if (node2 == node1){
      next
    }
    for (metric in metric_types){
      for (j in 1:length(tau)){
        metrics$metric_type[k] <- metric
        metrics$node_true[k] <- node1
        metrics$node_est[k] <- node2
        metrics$metric[k] <- 
          get_metric(spod_trends[ , paste(node2, tau[j], sep = "_")], 
                     spodPIDs[, node2], 
                     spod_signal[, node1], 
                     PID_thresh$thresh[PID_thresh$node == node2], 
                     metric)
        metrics$tau[k] <- tau[j]
        k <- k+1
        metrics[k+1, ] <- NA
      }
    }
  }
}

metric_choice <- "CAA"
metric_out <- metrics %>% 
  arrange(metric_type) %>%
  filter(tau == 0.1) %>%
  mutate(metric2 = round(metric,2))%>%
  select(node_true, node_est, metric2)

latex(metric_out, 
      file = "../Manuscript/application_metrics.tex",
      rowname = "",
      rowlabel = "",
      rgroup = metric_types,
      colheads = c("True", "Fitted"),
      n.rgroup = c(6, 6, 6, 6, 6))

spod_signal <- get_spod_signal(0.1, spod_trends, spodPIDs, PID_thresh$thresh)
spod_signal$time <- spodPIDs$time
spod_signal <- spod_signal %>%  gather(node, PID, -time) %>% filter(!is.na(node))
spodLong <- spodPeaks %>% gather("node","PID", -time) %>% filter(!is.na(node))

spodLong$node <- factor(spodLong$node)
spodLong_signal <- spodLong[which(spod_signal$PID==1), ]

ggplot(spodLong, aes(x=time, y=PID, col=node)) + 
  geom_line(alpha=0.5) + 
  geom_rug(data = spodLong_signal, sides = "b") +
  theme_bw() +
  facet_grid(node~.) + 
  xlim(c(as.POSIXct("2017-11-30 11:05:00"), 
         as.POSIXct("2017-11-30 11:30:00"))) + 
  ylim(c(-0.1, 5)) 
ggsave("../Manuscript/Figures/corrected_rugplot.png", width = 7, height = 3)

spod_signal[which(spod_signal$node == "f"), "PID"] <- 
  spod_signal[which(spod_signal$node == "f"), "PID"]*2
spod_signal[which(spod_signal$node == "g"), "PID"] <- 
  spod_signal[which(spod_signal$node == "g"), "PID"]*3

ggplot(spod_signal, aes(x=time, y=PID, col=node)) + 
  geom_line() + 
  theme_bw() +
  xlim(c(as.POSIXct("2017-11-30 11:05:00"), 
         as.POSIXct("2017-11-30 11:30:00"))) + 
  ylim(c(0, 3.5)) 
ggsave("../Manuscript/Figures/signal_plot.png", width = 7, height = 3)
