library(devtools)
load_all("detrendr")
library(tidyverse)

rm(list=ls())
# Application
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

plot(pid~time, spodNode[35000:37000, ], type="l")

y <- spodNode[1:1300, "pid"]
x <- spodNode[1:1300, "time"]
df.data <- data.frame(y=y, time=x)

n <- length(y)
k <- 3
tau <- c(0.05, 0.1)
overlap <- 100
rho <- 1
window_size <- 500
lambda <- 10*length(y)
max_iter <- 150
result0 <- gurobi_trend(y, tau, lambda, k)
result1 <- gurobi_trend(y[1:500], tau, lambda, k)
result2 <- gurobi_trend(y[401:900], tau, lambda, k)
result3 <- gurobi_trend(y[801:1300], tau, lambda, k)

df.sep.no <- data.frame(rbind(data.frame(time = x, result0, method="Single Fit"), 
                 data.frame(time = x[1:500], result1, method = "Window 1"),
                 data.frame(time = x[401:900], result2, method = "Window 2"),
                 data.frame(time = x[801:1300], result3, method = "Window 3"))) 
  

ggplot(df.data, aes(x=time, y=y)) +
  geom_line(col="grey") +
  geom_line(data = df.sep.no, aes(y=X1, col = "0.05", linetype = method))+
  geom_line(data = df.sep.no, aes(y=X2, col = "0.10", linetype = method))+
  scale_color_brewer(palette = "Set1")+
  labs(col="Quantile", linetype = "", x = "") + 
  theme_bw() 
ggsave("../Manuscript/overlapping_windows.png", width=7, height=1.5)

result <- consensus_ADMM(y, tau, lambda, k, rho, window_size, overlap, 
                         max_iter, eps = 2e-5, update = 1)


df_no <- rbind(data.frame(time=x, method = "Single Fit", result0), 
            data.frame(time=x, method = "Windows", result$theta))



ggplot(df.data, aes(x=time, y=y)) +
  geom_line(col="grey") +
  geom_line(data = df_no, aes(y=X1, col = "0.05", linetype = method))+
  geom_line(data = df_no, aes(y=X2, col = "0.10", linetype = method))+
  scale_color_brewer(palette = "Set1")+
  labs(col="Quantile", linetype = "", x = "") + 
  theme_bw() 
ggsave("../Manuscript/admm_windows.png", width=7, height=1.5)


y <- spodNode[34901:36200, "pid"]
x <- spodNode[34901:36200, "time"]
df.data <- data.frame(y=y, time=x)

n <- length(y)
k <- 3
tau <- c(0.05, 0.1)
overlap <- 100
rho <- 1
window_size <- 500
lambda <- 10*length(y)
max_iter <- 300
result0 <- gurobi_trend(y, tau, lambda, k)
result1 <- gurobi_trend(y[1:500], tau, lambda, k)
result2 <- gurobi_trend(y[401:900], tau, lambda, k)
result3 <- gurobi_trend(y[801:1300], tau, lambda, k)

df.sep <- data.frame(rbind(data.frame(time = x, result0, method="Single Fit"), 
                           data.frame(time = x[1:500], result1, method = "Window 1"),
                           data.frame(time = x[401:900], result2, method = "Window 2"),
                           data.frame(time = x[801:1300], result3, method = "Window 3"))) 


ggplot(df.data, aes(x=time, y=y)) +
  geom_line(col="grey") +
  geom_line(data = df.sep, aes(y=X1, col = "0.05", linetype = method))+
  geom_line(data = df.sep, aes(y=X2, col = "0.10", linetype = method))+
  scale_color_brewer(palette = "Set1")+
  labs(col="Quantile", linetype = "", x = "") + 
  theme_bw() 
ggsave("../Manuscript/overlapping_windows2.png", width=7, height=1.5)

result <- consensus_ADMM(y, tau, lambda, k, rho, window_size, overlap, 
                         max_iter, eps = 2e-5, update = 1)


df <- rbind(data.frame(time=x, method = "Single Fit", result0), 
            data.frame(time=x, method = "Windows", result$theta))



ggplot(df.data, aes(x=time, y=y)) +
  geom_line(col="grey") +
  geom_line(data = df, aes(y=X1, col = "0.05", linetype = method))+
  geom_line(data = df, aes(y=X2, col = "0.10", linetype = method))+
  scale_color_brewer(palette = "Set1")+
  labs(col="Quantile", linetype = "", x = "") + 
  theme_bw() 
ggsave("../Manuscript/admm_windows2.png", width=7, height=1.5)
