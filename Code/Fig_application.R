################################################################################
# Calculate mis-classification rates in peaks simulation
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
library(caret)
load_all("detrendr")
rm(list=ls())
spod <- read.csv("../SPod/fhrdata_2017-11-30.csv", 
                 header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

load("../SPod/spod_trends.RData")
nodes <- c("f", "g", "h")
spodPIDs <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")]/1000)
names(spodPIDs) <- nodes
spodPeaks <- spodPIDs - select(spod_trends, contains("0.1"))
spodPeaks$time <- spod_trends$time
spodPIDs$time <- spod_trends$time

spodRaw <- spod[, c("time", paste(nodes, "SPOD.PID..V.", sep="."))] %>%
  gather("node", "PID", -time)
ggplot(spodRaw, aes(x=time, y=PID, col=node)) +
  geom_line(alpha=0.5) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", labels = nodes)
ggsave("../Manuscript/Figures/uncorrected_data.png", width = 7, height = 3)

spodLong <- spodPeaks %>% gather("node","PID", -time)
thresholds <- 4*apply(spodPeaks[1:10000,1:3], 2, sd) +
  apply(spodPeaks[1:10000,1:3], 2, mean)
PID_thresh <- data.frame(node = names(thresholds), thresh = thresholds)

ggplot(spodLong, aes(x=time, y=PID, col=node)) +
  geom_line(alpha=0.5) +
  geom_hline(data=PID_thresh, aes(yintercept=thresh, col=node)) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", labels = nodes)
ggsave("../Manuscript/Figures/corrected_data.png", width = 7, height = 3)

# ################################################################################

get_spod_signal <- function(tau, spod_trends, spodPIDs, thresholds){
  spodPeaks <- select(spodPIDs, -time) - select(spod_trends, contains(paste(tau)))
  spodSignal <- spodPeaks
  for (i in 1:length(thresholds)){
    spodSignal[,i] <- as.numeric(spodPeaks[,i]>thresholds[i])
  }
  return(spodSignal)
}

spod_signal <- get_spod_signal(0.1, spod_trends, spodPIDs, PID_thresh$thresh)
spod_signal <- na.omit(spod_signal)
mat_h1 <- confusionMatrix(factor(spod_signal$f[spod_signal$h==1]),
                          factor(spod_signal$g[spod_signal$h==1]))

mat_h0 <- confusionMatrix(factor(spod_signal$f[spod_signal$h==0]),
                          factor(spod_signal$g[spod_signal$h==0]))

conf_out <- cbind(mat_h0$table, mat_h1$table)

latex(conf_out,
      file = "../Manuscript/complete_confusion.tex",
      rowlabel = "",
      rowname = c("f = 0", "f = 1"),
      cgroup = c("h = 0", "h = 1"),
      colheads = rep(c("g = 0", "g = 1"),2),
      n.cgroup = c(2,2),
      caption = "Confusion matrices for 3 SPod nodes after baseline removal.")

################################################################################
# Rug plot

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
