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
load("../SPod/spodPIDs.RData")
nodes <- c("c", "d", "e")
spodRaw <- spodPIDs %>%
  gather("node", "PID", -time)

ggplot(spodRaw, aes(x=time, y=PID, col=node)) +
  geom_line(alpha=0.5) +
  theme_bw() +
  #facet_grid(node~., scales = "free")+
  scale_color_brewer(palette = "Set1", labels = c("a", "b", "c"))+
  labs(col="SPod")
  # xlim(c(as.POSIXct("2017-04-13 12:30:01 EST"), 
  #        as.POSIXct("2017-04-13 13:30:01 EST")))
ggsave("../Manuscript/Figures/uncorrected_data.png", width = 7, height = 2.5)

load("../SPod/spod_trends.RData")
spodPeaks <- select(spodPIDs, -time) - select(spod_trends, contains("0.15"))
spodPeaks$time <- spod_trends$time
spodLong <- spodPeaks %>% gather("node","PID", -time)

ggplot(spodLong, aes(x=time, y=PID, col=node)) +
  geom_line(alpha=0.5) +
  theme_bw() +  
  scale_color_brewer(palette = "Set1", labels = c("a", "b", "c"))+
  labs(col="SPod")
ggsave("../Manuscript/Figures/corrected_data.png", width = 7, height = 2.5)

spodLong %>% filter(PID > 0.75)
# ################################################################################
detrendr_trends <- spod_trends
methods <- c("detrendr")
metric_df <- tibble(metric=NA, method=NA, tau=NA, crit=NA, metric_type=NA)
i <- 1
metrics <- c("confusion", "NMI", "VI")
nodes <- c("c", "d", "e")

for (method in methods){
  trends <- get(paste(method, "trends", sep = "_"))
  for (j in 1:length(tau)){
    for (crit in c(3, 4, 5)){
      signal <- get_spod_signal(tau[j], trends, spodPIDs, crit)
      for (metric in metrics){
        if (metric == "confusion"){
          metric_df$metric[i] <- I(list(get_confusion(signal, nodes)))
        } else if (metric == "NMI"){
          metric_df$metric[i] <- I(list(get_NMI(signal, nodes)))
        } else if (metric == "VI") {
          metric_df$metric[i] <- I(list(get_VI(signal, nodes)))
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
  filter(tau == 0.15, crit == 5, metric_type == "confusion") %>%
  select(metric) %>% 
  unnest() %>% gather("Label", "Value") %>%
  mutate(c = substr(Label, 6, 6), 
         d = substr(Label, 7, 7),
         e = substr(Label, 8, 8)) %>%
  arrange(c, d, e) %>% 
  select(-Label) %>%
  spread(c, Value) 

conf <- bind_cols(confusion_df[1:2,3:4], 
                          confusion_df[3:4,3:4])

metric_df %>% filter(tau == 0.15, crit == 5, metric_type == "VI") %>% unnest()

latex(conf,
      file = "../Manuscript/full_confusion_detrend.tex",
      rowlabel = "",
      rowname = c("b = 0", "b = 1", "b=0", "b=1"),
      cgroup = c("a = 0", "a = 1"),
      colheads = rep(c("c = 0", "c = 1"),2),
      n.cgroup = c(2,2),
      caption = "Confusion matrices for 3 SPods after baseline 
      removal using 15th quantile and threshold of 5*MAD (n=86,401).")

################################################################################
