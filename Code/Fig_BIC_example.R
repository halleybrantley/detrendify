library(tidyverse)
library(devtools)
load_all("detrendr")
rm(list=ls())

set.seed(9856882)
n <- 500
x <- seq(0.5, n, 1)/n
f <- sin(6*pi*x) 
y <- f + ((1+x^2)/4)*rnorm(n)
df <- data.frame(y=y, x=x, f=f)

y <- df$y
tau <- 0.1
k <- 3
lambdaSeq = length(y)^seq(0, 2, length.out=100)
df_trend <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
SIC <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
SIC2 <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
SIC3 <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
SIC4 <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
tau_min <- 0.5 - abs(0.5-tau)
D <- get_Dk(length(y), k)
gamma <- 1
df_tol <- 1e-9
for (i in 1:length(lambdaSeq)){
  lam <- lambdaSeq[i]
  suppressMessages(f_trend <- gurobi_trend(y, tau, lam, k))
  resid_trend <- checkloss(y-f_trend, tau)
  discr_diff <- abs(D%*%f_trend)
  df_trend[i,] <- Matrix::colSums(discr_diff > df_tol) #colSums(abs(y-f_trend)<df_tol)
  SIC[i] <- 2*colSums(resid_trend)/tau_min + log(n)*df_trend[i,] +
    2*gamma*log(choose(n-k, df_trend[i,]))
  SIC2[i] <- 2*colSums(resid_trend) + log(n)*df_trend[i,] 
  SIC3[i] <- 2*colSums(resid_trend)/tau_min + log(n)*df_trend[i,] 
  SIC4[i] <- log(colMeans(resid_trend)) + log(n)*df_trend[i,]/(2*n) 
}

plot.df <- data.frame(df = df_trend[,1])
plot.df$eBIC <- SIC
plot.df$BIC_noScale <- SIC2
plot.df$BIC_scale <- SIC3
plot.df$SIC <- n*(SIC4+5)
plot.df$lambda <- lambdaSeq



ggplot(plot.df, aes(x=log(lambda), y=df)) +
  geom_line() +
  theme_bw()
ggsave("../Manuscript/Figures/df_by_lambda.png", width = 3, height = 3)

plot.df %>% gather("criteria", "value", -c(df, lambda)) %>%
ggplot(aes(x=log(lambda), y=value, col = criteria)) +
  geom_line() + 
  theme_bw() +
  scale_color_brewer(palette = "Set1", 
                     labels = c("BIC no scale", "BIC with scale", "SIC",
                                "extended BIC")) +
  labs(y = "BIC", col="") +
  geom_vline(aes(xintercept = log(lambdaSeq[which.min(SIC)]), 
                 col = "eBIC")) +
  geom_vline(aes(xintercept = log(lambdaSeq[which.min(SIC3)]), 
                 col = "BIC_scale")) +
  geom_vline(aes(xintercept = log(lambdaSeq[which.min(SIC2)]), 
                 col = "BIC_noScale")) + 
  geom_vline(aes(xintercept = log(lambdaSeq[which.min(SIC4)]), 
                 col = "SIC")) 
ggsave("../Manuscript/Figures/BIC_by_lambda.png", width = 6, height = 3)


trend_eBIC <- gurobi_trend(y, tau, lambdaSeq[which.min(SIC)], k)
trend_BIC_no <- gurobi_trend(y, tau, lambdaSeq[which.min(SIC2)], k)
trend_BIC <- gurobi_trend(y, tau, lambdaSeq[which.min(SIC3)], k)
trend_SIC <- gurobi_trend(y, tau, lambdaSeq[which.min(SIC4)], k)
df$eBIC <- trend_eBIC
df$BIC_noScale <- trend_BIC_no
df$BIC_scale <- trend_BIC
df$SIC <- trend_SIC

ggplot(df, aes(x=x, y=y)) + 
  geom_line(col="grey") + 
  theme_bw() +
  geom_line(aes(y = SIC, col="SIC")) +
  geom_line(aes(y = BIC_scale, col = "BIC_scale")) +
  geom_line(aes(y = BIC_noScale, col = "BIC_noScale")) +
  geom_line(aes(y = eBIC, col="eBIC")) +
  scale_color_brewer(palette = "Set1", 
                     labels = c("BIC no scale", "BIC with scale", "SIC",
                                "extended BIC")) +
  labs(col = "") +
  guides(col = "none")
ggsave("../Manuscript/Figures/BIC_data.png", width = 3, height = 3)
