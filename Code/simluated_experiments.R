
library(Rglpk)
library(gurobi)
library(plyr)
library(devtools)
library(microbenchmark)
load_all("detrendr")

rm(list=ls())

n <- 400
x <- seq(1, n, 1)
y <- sin(x*2*pi/n) + rnorm(n, 0, .5)
plot(y, type="l")


k <- 3
tau <- c(0.05, 0.5)

lambdaSeq <- seq(50, 1600, 50)
loss_cv <- lambda_cv(y, tau, k, lambdaSeq,
                     numFolds = 5, parallel = TRUE)
plot(rowSums(loss_cv$loss))
loss_cv$lambda


theta_df <- get_windows(y, x, k, tau, loss_cv$lambda, n, 0)
plot.df <- left_join(theta_df, 
                     data.frame(time=x, y=y, 
                                true = sin(x*2*pi/n))) 

ggplot(plot.df, aes(x = time, y = y)) +
  geom_line(alpha = 0.2) +
  geom_line(aes(y=true), col="blue")+
  geom_line(aes(y=theta, col = factor(window), linetype=tau)) +
  theme_bw()


