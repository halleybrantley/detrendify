library(devtools)
library(doParallel)
build("detrendr")
install("detrendr")
library(detrendr)
library(tidyverse)

rm(list=ls())
# Application
dataDir <- "~/Desktop/EPA/SPod_Data/PHL_2017"
datafiles <- dir(dataDir, pattern=".csv", full.names=TRUE)
spod <- read.csv(datafiles[2], header=TRUE,  na.strings = "N/A")
spod$time <- as.POSIXct(strptime(as.character(spod$TimeStamp), 
                                 format= "%m/%d/%Y %H:%M:%S")) 


node <- "c"
pidCol <- paste(node, "SPOD.PID..V.", sep=".")
spodNode <- spod[, c("time", pidCol)]
names(spodNode)[2] <- c("pid")
spodNode <- subset(spodNode, !is.na(pid))
spodNode$pid <- as.numeric(scale(spodNode$pid, center = TRUE))
plot(pid~time, spodNode, type="l")

y <- spodNode[, "pid"]
x <- spodNode[, "time"]
k <- 3
tau <- c(0.1, 0.5)
overlap <- 300
window_size <- 3600
lambda <- 10*window_size

theta.df <- get_windows(y, x, k, tau, lambda, window_size, overlap)
plot.df <- left_join(theta.df, spodNode) 

ggplot(plot.df, aes(x = time, y = pid)) +
  geom_line(alpha = 0.2) +
  geom_line(aes(y=theta, col = factor(window), linetype=tau)) +
  theme_bw() + 
  guides(col = "none") + 
  ggtitle("Finished in 64 s")

ggsave("windows_fit.png", width = 8, height = 4)








k <- 3
tau <- c(0.05, 0.5)
overlap <- 300
rho <- 1
window_size <- 2400
lambda <- 5*window_size
max_iter <- 5
cl <- makeCluster(detectCores() - 1 , type="FORK")
clusterEvalQ(cl, library(detrendr))

system.time(result <- consensus_ADMM(y, tau, lambda, k, rho, window_size, overlap, 
                                     max_iter, eps = 5e-4, update = 1, cl = cl))
# Vector of 40,000 window of 2400 overlap 300, converged in 3 interations and 1045 s

closeAllConnections()
showConnections()

y_n <- length(y)

plot(y, type="l", col="grey", main = "Converged in 3 interations and 1045 s")
lines(result$theta[,1], col="red")
lines(result$theta[,2], col="blue")
abline(v=window_size, col="darkgreen")
abline(v=window_size-overlap, col="purple")

n_windows <- ceiling(length(y)/(window_size-overlap))

for (i in 2:n_windows){
  abline(v=i*window_size-(i-1)*overlap, col="purple")
  abline(v=i*window_size-i*overlap, col="darkgreen")
}

dev.copy(png, "consensus_fig.png", width = 800, height = 400)
dev.off()

system.time(theta.df <- get_windows(y, x, k, tau, lambda, window_size, overlap))
plot.df <- left_join(theta.df, spodNode) 

ggplot(plot.df, aes(x = time, y = pid)) +
  geom_line(alpha = 0.2) +
  geom_line(aes(y=theta, col = factor(window), linetype=tau)) +
  theme_bw() + 
  guides(col = "none") + 
  ggtitle("Finished in 64 s")

ggsave("windows_fit.png", width = 8, height = 4)


y <- spodNode[12000:13000, "pid"]
x <- spodNode[12000:13000, "time"]
k <- 3
tau <- c(0.05, 0.5)
lambda <- 5*window_size
rho <- 1
w <- numeric(length(y)*length(tau))
phiBar <- numeric(length(y)*length(tau))
D <- as.matrix(get_Dk(length(y), k))
library(microbenchmark)

microbenchmark(
  theta1 <-  quad_update(y, tau, lambda, D, w, phiBar, rho),
  theta2 <- gurobi_trend(y, tau, lambda, k),
  times = 10
)
