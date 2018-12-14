library(devtools)
library(gurobi)
load_all("detrendr")
rm(list=ls())

dataDir <- "~/Desktop/EPA/SPod_Data/TestRange_Dec2017"
datafiles <- dir(dataDir, pattern=".csv", full.names=TRUE)

spod <- read.csv(datafiles[3], header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 

node <- "h"
pidCol <- paste(node, "SPOD.PID..V.", sep=".")
spodNode <- spod[, c("time", pidCol)]
names(spodNode)[2] <- c("pid")
spodNode <- subset(spodNode, !is.na(pid))
spodNode$pid <- as.numeric(scale(spodNode$pid, center = TRUE))
window_size <- 5000
overlap <- 500
max_iter <- 3
lambda <- 2*window_size
tau <- 0.1
k <- 3
spod_trends <- data.frame(time = spod$time, f=NA, g = NA, h = NA)

for (node in c("f", "g", "h")){
  pidCol <- paste(node, "SPOD.PID..V.", sep=".")
  spodNode <- spod[, c("time", pidCol)]
  names(spodNode)[2] <- c("pid")
  spodNode$pid <- as.numeric(scale(spodNode$pid, center = TRUE))
  result <- get_windows_BIC(spodNode$pid, tau, k, window_size, overlap,
                          lambdaSeq = window_size^seq(1.2, 1.5, length.out=5), 
                          df_tol = 1e-9, 
                          gamma = 1,
                          plot_lambda = FALSE, 
                          solver = NULL, 
                          criteria = "eBIC")
  spod_trends[,node] <- result$trend
}

nodes <- c("f", "g", "h")
spodPIDs <- as.data.frame(scale(spod[, paste(nodes, "SPOD.PID..V.", sep=".")], 
                                center=TRUE))
names(spodPIDs) <- nodes
spodPeaks <- spodPIDs - spod_trends[,2:4]

spodPeaks$time <- spod_trends$time  
spodPIDs$time <- spod_trends$time

spodRaw <- spodPIDs %>% gather("node", "PID", -time)
ggplot(spodRaw, aes(x=time, y=PID, col=node)) + 
  geom_line(alpha=0.5) + 
  theme_bw() 
ggsave("../Manuscript/Figures/uncorrected_data.png", width = 7, height = 3)

spodLong <- spodPeaks %>% gather("node","PID", -time)

ggplot(spodLong, aes(x=time, y=PID, col=node)) + 
  geom_line(alpha=0.5) + 
  theme_bw() 
ggsave("../Manuscript/Figures/corrected_data.png", width = 7, height = 3)

ggplot(spodLong, aes(x=time, y=PID, col=node)) + 
  geom_line(alpha=0.5) + 
  theme_bw() +
  xlim(c(as.POSIXct("2017-11-30 9:09:57"), 
         as.POSIXct("2017-11-30 12:29:57")))
ggsave("../Manuscript/Figures/corrected_zoom_data.png", width = 7, height = 3)

plot(spodNode$pid, type="l")
lines(result$trend[,1], col="red")

trend <-  get_trend_windows(spodNode$pid, tau, lambda, k, rho=1, window_size,
                             overlap, max_iter, update=1, 
                             quad = TRUE)
plot(trend2, type="l")
lines(trend, type="l", col="red")
