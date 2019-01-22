################################################################################
# Calculate mis-classification rates in peaks simulation
# Halley Brantley
################################################################################
library(tidyverse)
library(devtools)
library(Hmisc)
library(caret)
load_all("detrendr")
source("application_functions.R")
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

spodRaw <- spodPIDs %>%
  gather("node", "PID", -time)
ggplot(spodRaw, aes(x=time, y=PID, col=node)) +
  geom_line(alpha=0.5) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", labels = nodes)
ggsave("../Manuscript/Figures/uncorrected_data.png", width = 7, height = 2.5)

spodLong <- spodPeaks %>% gather("node","PID", -time)

ggplot(spodLong, aes(x=time, y=PID, col=node)) +
  geom_line(alpha=0.5) +
  theme_bw() +
  scale_color_brewer(palette = "Set1", labels = nodes)
ggsave("../Manuscript/Figures/corrected_data.png", width = 7, height = 2.5)

# ################################################################################

NMI <- {}
conf <- {}
for (j in 1:length(tau)){
  detrendr_signal <- get_spod_signal(tau[j], spod_trends, spodPIDs, crit = 7)

  conf <- rbind(conf, get_confusion(detrendr_signal))
  NMI <- rbind(NMI, 
               cbind("fg", round(NMI(detrendr_signal$f, detrendr_signal$g, 
                                     variant="sqrt"),2)), 
               cbind("fh", round(NMI(detrendr_signal$f, detrendr_signal$h, 
                                     variant="sqrt"),2)), 
               cbind("gh", round(NMI(detrendr_signal$g, detrendr_signal$h, 
                                     variant="sqrt"),2)) 
  )
}

  
latex(conf,
      file = "../Manuscript/full_confusion_detrend.tex",
      rowlabel = "",
      rowname = c("f = 0", "f = 1", "f=0", "f=1"),
      cgroup = c("h = 0", "h = 1"),
      rgroup = c("tau=0.05", "tau=0.1"),
      n.rgroup = c(2, 2),
      colheads = rep(c("g = 0", "g = 1"),2),
      n.cgroup = c(2,2),
      caption = "Confusion matrices for 3 SPod nodes after baseline removal (n=52322).")

latex(NMI[,2], 
      file = "../Manuscript/full_NMI.tex", 
      rowlabel="",
      rowname = NMI[,1],
      n.rgroup = c(3,3), 
      rgroup = c("tau=0.05", "tau=0.1"), 
      caption = "Normalized mutual information between signal classifications, full dataset."
)

################################################################################
signal <- get_spod_signal(tau[1], spod_trends, spodPIDs, crit = 7)
spodPeaks <- select(spodPIDs, -time) - select(spod_trends, contains("0.05"))
spodPeaks$time <- spodPIDs$time
max(which(apply(signal==1, 1, any)))

signal <- signal[34073:49193, ]
spod_sub <- spodPeaks[34073:49193, ]
signal_long <- signal %>%  gather(node, PID) %>% filter(!is.na(node))
spodLong <- spod_sub %>% gather("node","PID", -time) %>% filter(!is.na(node))

spodLong$node <- factor(spodLong$node)
spodLong_signal <- spodLong[which(signal_long$PID==1), ]

spodPeaks[is.na(spodPeaks)] <- 0
crit <- 7
thresholds <- data.frame(node = nodes, 
                         threshold = apply(spodPeaks[,nodes], 2, function(x) 
                                             median(x, na.rm=T) + 
                                             crit*mean(abs(x-median(x, na.rm=T)), na.rm=T)), 
                         median = apply(spodPeaks[,nodes], 2, function(x) 
                           median(x, na.rm=T)), 
                         mad = apply(spodPeaks[,nodes], 2, function(x) 
                           median(x, na.rm=T) + 
                             median(abs(x-median(x, na.rm=T)), na.rm=T)))
                    
mad <- median(abs(spod_sub$f - median(spod_sub$f, na.rm=T)), na.rm=T)
med <- median(spod_sub$f, na.rm=T)
hist(spod_sub$f, 50)
abline(v=median(spod_sub$f, na.rm=T), col="red")
abline(v=median(spod_sub$f, na.rm=T)+mad, col="blue")
plot(spod_sub$f, type="l")
abline(h=med, col="red")
abline(h=med+mad, col="blue")

ggplot(spodLong, aes(x=time, y=PID)) +
  geom_rug(data = spodLong_signal, sides = "b") +
  geom_line(col="grey") +
  geom_hline(data = thresholds, aes(yintercept=threshold))+
  geom_hline(data = thresholds, aes(yintercept=median))+
  geom_hline(data = thresholds, aes(yintercept=mad))+
  theme_bw() +
  facet_grid(node~., scales="free_y")

