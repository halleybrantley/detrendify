library(fields)
library(tidyverse)
library(devtools)
library(Hmisc)
library(caret)
library(zoo)
library(aricode)
load_all("detrendr")
source("application_functions.R")
rm(list=ls())
load("../SPod/trends_short.RData")

spodPeaks <- select(spodPIDs, -time) - select(detrendr_trends,
                                              contains(paste(0.15)))
plot(spodPeaks$h, type="l")
thresholds <- apply(spodPeaks, 2, 
                      function(x) median(x, na.rm=T) + 
                      4*median(abs(x-median(x, na.rm=T)), na.rm=T))
abline(h=thresholds[3], col="red")
plot(g~f, spodPeaks)


################################################################################
crit <- 5
NMI <- {}
conf <- {}
for (j in 1:length(tau)){
  detrendr_signal <- get_spod_signal(tau[j], detrendr_trends, spodPIDs, crit)
  qsreg_signal <- get_spod_signal(tau[j], qsreg_trends, spodPIDs, crit)
  conf <- rbind(conf, cbind(get_confusion(detrendr_signal), 
                            get_confusion(qsreg_signal)))
  NMI <- rbind(NMI, 
               cbind("fg", round(NMI(detrendr_signal$f, detrendr_signal$g, 
                               variant="sqrt"),2),
                     round(NMI(qsreg_signal$f, qsreg_signal$g, 
                               variant="sqrt"), 2)), 
               cbind("fh", round(NMI(detrendr_signal$f, detrendr_signal$h, 
                               variant="sqrt"),2),
                     round(NMI(qsreg_signal$f, qsreg_signal$h, 
                               variant="sqrt"), 2)), 
               cbind("gh", round(NMI(detrendr_signal$g, detrendr_signal$h, 
                               variant="sqrt"),2),
                     round(NMI(qsreg_signal$g, qsreg_signal$h, 
                               variant="sqrt"), 2)) 
  )
}

conf <- rbind(rep(c("g = 0", "g = 1"),4), conf)

latex(conf,
      file = "../Manuscript/short_confusion_detrend.tex",
      rowlabel = "",
      rowname = c("", "f = 0", "f = 1", "f=0", "f=1"),
      cgroup = c("detrendr", "qsreg"),
      rgroup = c("", "tau=0.05", "tau=0.1"),
      n.rgroup = c(1, 2, 2),
      colheads = rep(c("h = 0", "", "h = 1", ""),2),
      n.cgroup = c(4,4),
      caption = "Confusion matrices for 3 SPod nodes after baseline removal (n=8000).")

latex(NMI[,2:3], 
      file = "../Manuscript/NMI.tex", 
      rowlabel="",
      rowname = NMI[,1],
      colheads = c("detrendr", "qsreg"), 
      n.rgroup = c(3,3), 
      rgroup = c("tau=0.05", "tau=0.1"), 
      caption = "Normalized mutual information between signal classifications."
)
################################################################################
# Rug plot

spod_signal <- get_spod_signal(0.1, detrendr_trends, spodPIDs)
spod_signal$time <- spodPIDs$time
spod_signal <- spod_signal %>%  gather(node, PID, -time) %>% filter(!is.na(node))
spodPeaks$time <- spodPIDs$time
spodLong <- spodPeaks %>% gather("node","PID", -time) %>% filter(!is.na(node))

spodLong$node <- factor(spodLong$node)
spodLong_signal <- spodLong[which(spod_signal$PID==1), ]

ggplot(spodLong, aes(x=time, y=PID)) +
  geom_rug(data = spodLong_signal, sides = "b") +
  geom_line(col="grey") +
  theme_bw() +
  facet_grid(node~., scales="free_y")
ggsave("../Manuscript/Figures/corrected_rugplot.png", width = 7, height = 3)
