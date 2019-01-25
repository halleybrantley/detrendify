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
source("application_functions.R")
spod <- read.csv("../SPod/fhrdata_2017-11-30.csv", 
                 header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

load("../SPod/spod_trends.RData")

nodes <- c("f", "g", "h")
spodPIDs <- as.data.frame(spod[, paste(nodes, "SPOD.PID..V.", sep=".")]/1000)
names(spodPIDs) <- nodes
spodPeaks <- spodPIDs - select(spod_trends, contains("0.15"))
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
detrendr_trends <- spod_trends
methods <- c("detrendr")
metric_df <- tibble(metric=NA, method=NA, tau=NA, crit=NA, metric_type=NA)
i <- 1
metrics <- c("confusion", "NMI", "VI")

for (method in methods){
  trends <- get(paste(method, "trends", sep = "_"))
  for (j in 1:length(tau)){
    for (crit in c(3, 4, 5)){
      signal <- get_spod_signal(tau[j], trends, spodPIDs, crit)
      for (metric in metrics){
        if (metric == "confusion"){
          metric_df$metric[i] <- I(list(get_confusion(signal)))
        } else if (metric == "NMI"){
          metric_df$metric[i] <- I(list(get_NMI(signal)))
        } else if (metric == "VI") {
          metric_df$metric[i] <- I(list(get_VI(signal)))
        }
        metric_df$method[i] <- method
        metric_df$tau[i] <- tau[j]
        metric_df$crit[i] <- crit
        metric_df$metric_type[i] <- metric
        i <- i+1
        metric_df[i, ] <- NA
      }
    }
  }
}
metric_df <- metric_df[-i,]

confusion_df <- metric_df %>% 
  filter(tau == 0.15, crit == 3, metric_type == "confusion") %>%
  select(metric) %>% 
  unnest() %>% gather("Label", "Value") %>%
  mutate(h = substr(Label, 2, 2), 
         f = substr(Label, 4, 4),
         g = substr(Label, 6, 6)) %>%
  arrange(f, g, h) %>% 
  select(-Label) %>%
  spread(h, Value) 

conf <- bind_cols(confusion_df[1:2,3:4], 
                          confusion_df[3:4,3:4])

metric_df %>% filter(tau == 0.15, crit == 3, metric_type == "VI") %>% unnest()

latex(conf,
      file = "../Manuscript/full_confusion_detrend.tex",
      rowlabel = "",
      rowname = c("g = 0", "g = 1", "g=0", "g=1"),
      cgroup = c("f = 0", "f = 1"),
      colheads = rep(c("h = 0", "h = 1"),2),
      n.cgroup = c(2,2),
      caption = "Confusion matrices for 3 SPod nodes after baseline 
      removal using 15th quantil and threshold of 3*MAD (n=52322).")

################################################################################
