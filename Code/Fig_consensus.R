library(tidyverse)
library(devtools)
load_all("detrendr")

rm(list=ls())
source("sim_generating_functions.R")
set.seed(987652)
overlap <- 150
window_size <- 500
n <- window_size*3 - overlap*2
df.data <- generate_peaks_design(n)
df.data$x <- seq(1, n, 1)
y <- df.data$y
x <- df.data$x
k <- 3
tau <- c(0.05, 0.1)

lambda <- length(y)
max_iter <- 25
result0 <- get_trend(y, tau, lambda, k)
result1 <- get_trend(y[1:window_size], tau, lambda, k)
result2 <- get_trend(y[(window_size - overlap + 1):(2*window_size - overlap)], tau, lambda, k)
result3 <- get_trend(y[(2*window_size - 2*overlap + 1):(3*window_size - 2*overlap)], tau, lambda, k)

df.sep.no <- data.frame(rbind(data.frame(x = x, result0, method="Single Fit"), 
                              data.frame(x = x[1:window_size], result1, method = "Window 1"),
                              data.frame(x = x[(window_size - overlap + 1):(2*window_size - overlap)], result2, method = "Window 2"),
                              data.frame(x = x[(2*window_size - 2*overlap + 1):(3*window_size - 2*overlap)], result3, method = "Window 3"))) 


ggplot(df.data, aes(x=x, y=y)) +
  geom_line(col="grey") +
  geom_line(data = df.sep.no, aes(y=X1, col = "0.05", linetype = method))+
  geom_line(data = df.sep.no, aes(y=X2, col = "0.10", linetype = method))+
  scale_color_brewer(palette = "Set1")+
  labs(col="Quantile", linetype = "", x = "") + 
  theme_bw() +
  geom_segment(aes(x = x[500], xend=x[500], y=0, yend=1.5)) +
  geom_label(aes(x=x[500], y=-.2, label = bquote(u_1)))
ggsave("../Manuscript/Figures/overlapping_windows.png", width=7, height=2.5)

result <- get_trend_windows(y, tau, lambda, k, window_size, overlap, 
                            max_iter, update = 1, eps_abs = .01)


df_no <- rbind(data.frame(x=x, method = "Single Fit", result0), 
               data.frame(x=x, method = "Windows", result))


ggplot(df.data, aes(x=x, y=y)) +
  geom_line(col="grey") +
  geom_line(data = df_no, aes(y=X1, col = "0.05", linetype = method))+
  geom_line(data = df_no, aes(y=X2, col = "0.10", linetype = method))+
  scale_color_brewer(palette = "Set1")+
  labs(col="Quantile", linetype = "", x = "") + 
  theme_bw() 
ggsave("../Manuscript/Figures/admm_windows.png", width=7, height=2.5)
max(abs(result0-result))

################################################################################
# Application
# dataDir <- "~/Desktop/EPA/SPod_Data/TestRange_Dec2017"
# datafiles <- dir(dataDir, pattern=".csv", full.names=TRUE)
# spod <- read.csv(datafiles[3], header=TRUE,  na.strings = "N/A")
# spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
#                                  format= "%m/%d/%Y %H:%M:%S"))
# 
# node <- "h"
# pidCol <- paste(node, "SPOD.PID..V.", sep=".")
# spodNode <- spod[, c("time", pidCol)]
# names(spodNode)[2] <- c("pid")
# spodNode <- subset(spodNode, !is.na(pid))
# spodNode$pid <- as.numeric(scale(spodNode$pid, center = TRUE))
# 
# plot(pid~time, spodNode[35000:37000, ], type="l")
# 
# overlap <- 200
# window_size <- 500
# y <- spodNode[1:(3*window_size - 2*overlap), "pid"]
# x <- spodNode[1:(3*window_size - 2*overlap), "time"]
# df.data <- data.frame(y=y, time=x)
# 
# n <- length(y)
# k <- 3
# tau <- c(0.2, 0.3)
# rho <- 3
# lambda <- 10*length(y)
# max_iter <- 25
# result0 <- get_trend(y, tau, lambda, k)
# result1 <- get_trend(y[1:window_size], tau, lambda, k)
# result2 <- get_trend(y[(window_size - overlap + 1):(2*window_size - overlap)], tau, lambda, k)
# result3 <- get_trend(y[(2*window_size - 2*overlap + 1):(3*window_size - 2*overlap)], tau, lambda, k)
# 
# df.sep.no <- data.frame(rbind(data.frame(time = x, result0, method="Single Fit"), 
#                  data.frame(time = x[1:window_size], result1, method = "Window 1"),
#                  data.frame(time = x[(window_size - overlap + 1):(2*window_size - overlap)], result2, method = "Window 2"),
#                  data.frame(time = x[(2*window_size - 2*overlap + 1):(3*window_size - 2*overlap)], result3, method = "Window 3"))) 
#   
# 
# ggplot(df.data, aes(x=time, y=y)) +
#   geom_line(col="grey") +
#   geom_line(data = df.sep.no, aes(y=X1, col = "0.05", linetype = method))+
#   geom_line(data = df.sep.no, aes(y=X2, col = "0.10", linetype = method))+
#   scale_color_brewer(palette = "Set1")+
#   labs(col="Quantile", linetype = "", x = "") + 
#   theme_bw() 
# ggsave("../Manuscript/Figures/app_overlapping_windows.png", width=7, height=1.5)
# 
# result <- get_trend_windows(y, tau, lambda, k, rho, window_size, overlap, 
#                          max_iter, update = 5)
# 
# 
# df_no <- rbind(data.frame(time=x, method = "Single Fit", result0), 
#             data.frame(time=x, method = "Windows", result))
# 
# 
# 
# ggplot(df.data, aes(x=time, y=y)) +
#   geom_line(col="grey") +
#   geom_line(data = df_no, aes(y=X1, col = "0.05", linetype = method))+
#   geom_line(data = df_no, aes(y=X2, col = "0.10", linetype = method))+
#   scale_color_brewer(palette = "Set1")+
#   labs(col="Quantile", linetype = "", x = "") + 
#   theme_bw() 
# ggsave("../Manuscript/Figures/app_admm_windows.png", width=7, height=1.5)
# 
# 
# y <- spodNode[34901:36000, "pid"]
# x <- spodNode[34901:36000, "time"]
# df.data <- data.frame(y=y, time=x)
# overlap <- 200
# n <- length(y)
# k <- 3
# tau <- c(0.05, 0.1)
# rho <- 1
# lambda <- 10*length(y)
# max_iter <- 25
# result0 <- get_trend(y, tau, lambda, k)
# result1 <- get_trend(y[1:window_size], tau, lambda, k)
# result2 <- get_trend(y[(window_size - overlap + 1):(2*window_size - overlap)], tau, lambda, k)
# result3 <- get_trend(y[(2*window_size - 2*overlap + 1):(3*window_size - 2*overlap)], tau, lambda, k)
# 
# df.sep <- data.frame(rbind(data.frame(time = x, result0, method="Single Fit"), 
#                            data.frame(time = x[1:window_size], result1, method = "Window 1"),
#                            data.frame(time = x[(window_size - overlap + 1):(2*window_size - overlap)], result2, method = "Window 2"),
#                            data.frame(time = x[(2*window_size - 2*overlap + 1):(3*window_size - 2*overlap)], result3, method = "Window 3"))) 
# 
# 
# ggplot(df.data, aes(x=time, y=y)) +
#   geom_line(col="grey") +
#   geom_line(data = df.sep, aes(y=X1, col = "0.05", linetype = method))+
#   geom_line(data = df.sep, aes(y=X2, col = "0.10", linetype = method))+
#   scale_color_brewer(palette = "Set1")+
#   labs(col="Quantile", linetype = "", x = "") + 
#   theme_bw() 
# ggsave("../Manuscript/Figures/overlapping_windows2.png", width=7, height=2.5)
# 
# result <- get_trend_windows(y, tau, lambda, k, rho, window_size, overlap, 
#                          max_iter, update = 1)
# 
# 
# df <- rbind(data.frame(time=x, method = "Single Fit", result0), 
#             data.frame(time=x, method = "Windows", result))
# 
# 
# 
# ggplot(df.data, aes(x=time, y=y)) +
#   geom_line(col="grey") +
#   geom_line(data = df, aes(y=X1, col = "0.05", linetype = method))+
#   geom_line(data = df, aes(y=X2, col = "0.10", linetype = method))+
#   scale_color_brewer(palette = "Set1")+
#   labs(col="Quantile", linetype = "", x = "") + 
#   theme_bw() 
# ggsave("../Manuscript/Figures/admm_windows2.png", width=7, height=2.5)
