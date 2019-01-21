library(fields)
library(tidyverse)
library(devtools)
library(Hmisc)
library(caret)
library(zoo)
library(aricode)
load_all("detrendr")
rm(list=ls())
load("../SPod/trends_short.RData")

spodPeaks <- select(spodPIDs, -time) - select(detrendr_trends,
                                              contains(paste(0.05)))
plot(spodPeaks$f, type="l")
thresholds <- apply(spodPeaks, 2, 
                      function(x) median(x) + 2.5*mean(abs(x-median(x))))
abline(h=thresholds[1], col="red")


################################################################################
get_spod_signal <- function(tau, spod_trends, spodPIDs){
  spodPeaks <- select(spodPIDs, -time) - select(spod_trends, contains(paste(tau)))
  thresholds <- apply(spodPeaks, 2, 
                  function(x) median(x) + 3*mean(abs(x-median(x))))
  spodSignal <- spodPeaks
  for (i in 1:length(thresholds)){
    spodSignal[,i] <- as.numeric(spodPeaks[,i]>thresholds[i])
  }
  return(spodSignal)
}

get_confusion <- function(spod_signal){
  spod_signal <- na.omit(spod_signal)
  for (i in 1:ncol(spod_signal)){
    spod_signal[,i] <- factor(spod_signal[,i])
  }
  mat_h1 <- confusionMatrix(spod_signal$f[spod_signal$h==1],
                            spod_signal$g[spod_signal$h==1])

  mat_h0 <- confusionMatrix(spod_signal$f[spod_signal$h==0],
                            spod_signal$g[spod_signal$h==0])

  conf_out <- cbind(mat_h0$table, mat_h1$table)
  return(conf_out)
}

NMI <- {}
conf <- {}
for (j in 1:length(tau)){
  detrendr_signal <- get_spod_signal(tau[j], detrendr_trends, spodPIDs)
  qsreg_signal <- get_spod_signal(tau[j], qsreg_trends, spodPIDs)
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
