################################################################################
# Extended and Scaled BIC Illustration Figures
################################################################################
library(tidyverse)
library(devtools)
load_all("../../detrendr")
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
  f_trend <- get_trend(y, tau, lam, k)
  resid_trend <- checkloss(y-f_trend, tau)
  discr_diff <- abs(D%*%f_trend)
  df_trend[i,] <- Matrix::colSums(discr_diff > df_tol) 
  SIC[i] <- 2*colSums(resid_trend)/tau_min + log(n)*df_trend[i,] +
    2*gamma*log(choose(n-k, df_trend[i,]))
  SIC2[i] <- 2*colSums(resid_trend) + log(n)*df_trend[i,] 
  SIC3[i] <- 2*colSums(resid_trend)/tau_min + log(n)*df_trend[i,] 
  SIC4[i] <- log(colMeans(resid_trend)) + log(n)*df_trend[i,]/(2*n) 
}

plot.df <- data.frame(df = df_trend[,1])
plot.df$eBIC <- SIC-min(SIC)
plot.df$BIC_noScale <- SIC2-min(SIC2)
plot.df$BIC_scale <- SIC3 - min(SIC3)
plot.df$SIC <- n*(SIC4 - min(SIC4))
plot.df$lambda <- lambdaSeq



ggplot(plot.df, aes(x=log(lambda), y=df)) +
  geom_line() +
  theme_bw() +
  theme(text = element_text(size=22))+
  labs(x=bquote(log(lambda)))
  ggsave("../Manuscript/Figures/df_by_lambda.png", width = 3, height = 3)

colPal <- c('#006d2c', '#2ca25f', '#66c2a4', '#1c9099', '#d8b365',
              "#c2a5cf", "#9970ab", "#762a83")

BIC.df <- plot.df %>% gather("criteria", "value", -c(df, lambda)) %>%
  mutate(criteria = factor(criteria, 
                           levels = c("BIC_noScale", "BIC_scale", "eBIC", "SIC")))

ggplot(BIC.df, aes(x=log(lambda), y=value, col = criteria)) +
  geom_line(size=1.2) + 
  theme_bw() +
  scale_color_manual(values = c('#fd8d3c', '#df65b0', '#006d2c', '#2ca25f'), 
                     labels = c("BIC no scale", "BIC with scale",
                                "eBIC with scale", "SIC")) +
  labs(y = "BIC (shifted)", col="") +
  geom_vline(aes(xintercept = log(lambdaSeq[which.min(SIC)]), 
                 col = "eBIC"), linetype=2) +
  geom_vline(aes(xintercept = log(lambdaSeq[which.min(SIC3)]), 
                 col = "BIC_scale"), linetype=2) +
  geom_vline(aes(xintercept = log(lambdaSeq[which.min(SIC2)]), 
                 col = "BIC_noScale"), linetype=2) + 
  geom_vline(aes(xintercept = log(lambdaSeq[which.min(SIC4)]), 
                 col = "SIC"), linetype=2) +
  labs(x=bquote(log(lambda))) +
  theme(text = element_text(size=18))
ggsave("../Manuscript/Figures/BIC_by_lambda.png", width = 6, height = 4)


trend_eBIC <- get_trend(y, tau, lambdaSeq[which.min(SIC)], k)
trend_BIC_no <- get_trend(y, tau, lambdaSeq[which.min(SIC2)], k)
trend_BIC <- get_trend(y, tau, lambdaSeq[which.min(SIC3)], k)
trend_SIC <- get_trend(y, tau, lambdaSeq[which.min(SIC4)], k)
df$eBIC <- trend_eBIC
df$BIC_noScale <- trend_BIC_no
df$BIC_scale <- trend_BIC
df$SIC <- trend_SIC

ggplot(df, aes(x=x*100, y=y)) + 
  geom_line(col="grey") + 
  theme_bw() +
  geom_line(aes(y = SIC, col="SIC")) +
  geom_line(aes(y = BIC_scale, col = "BIC_scale")) +
  geom_line(aes(y = BIC_noScale, col = "BIC_noScale")) +
  geom_line(aes(y = eBIC, col="eBIC"), size=1.2) +
  theme(text = element_text(size=22)) +
  scale_color_manual(values = c("BIC_noScale"='#fd8d3c', 
                                   "BIC_scale"='#df65b0', 
                                   "eBIC"='#006d2c', "SIC"='#2ca25f'))  +
  labs(col = "", x = "x") +
  guides(col = "none")
ggsave("../Manuscript/Figures/BIC_data.png", width = 4, height = 4)
