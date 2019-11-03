################################################################################
# Extended and Scaled BIC Illustration Figures
################################################################################
library(tidyverse)
library(devtools)
load_all("../../detrendr")
rm(list=ls())

set.seed(2019)
n <- 500
x <- seq(0.5, n, 1)/n
f <- sin(6*pi*x) 
y <- f + ((1+x^2)/4)*rnorm(n)
df <- data.frame(y=y, x=x, f=f)

y <- df$y
tau <- 0.1
k <- 3
lambdaSeq = exp(seq(0, 8, length.out=100))
df_trend <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
sigs <- seq(.1, 1, .3)
SIC <- matrix(NA, nrow=length(lambdaSeq), ncol=length(sigs))
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
  for (j in 1:length(sigs)){
    SIC[i, j] <- 2*colSums(resid_trend)/sigs[j] + log(n)*df_trend[i,] +
      2*gamma*log(choose(n-k, df_trend[i,]))
  }
}

plot.df <- data.frame(df = df_trend[,1])
plot.df[,2:5] <- sweep(SIC,2,apply(SIC, 2, min), "-")
plot.df$lambda <- lambdaSeq
names(plot.df)[2:5] <- paste("sigma", sigs)
plot.df <- plot.df %>% filter(log(lambda) < 9.5)

ggplot(plot.df, aes(x=log(lambda), y=df)) +
  geom_line() +
  theme_bw() +
  theme(text = element_text(size=22))+
  labs(x=bquote(log(lambda)))

colPal <- c('#006d2c', '#2ca25f', '#66c2a4', '#1c9099', '#d8b365',
              "#c2a5cf", "#9970ab", "#762a83")

BIC.df <- plot.df %>% 
  gather("criteria", "value", -c(df, lambda)) 

ggplot(BIC.df, aes(x=log(lambda), y=value, col = criteria)) +
  geom_line() + 
  theme_bw() +
  scale_color_manual(values = c('#fd8d3c', '#df65b0', '#006d2c', '#2ca25f')) +
  labs(y = "eBIC (shifted)", col="")
ggsave("../Manuscript/Figures/BIC_by_lambda_and_sigma.png", width = 6, height = 4)


trend_eBIC <- get_trend(y, tau, lambdaSeq[which.min(SIC[,1])], k)
df$eBIC <- trend_eBIC

ggplot(df, aes(x=x*100, y=y)) + 
  geom_line(col="grey") + 
  theme_bw() +
  geom_line(aes(y = eBIC, col="eBIC"), size=1.2) +
  theme(text = element_text(size=22)) +
  labs(col = "", x = "x") +
  guides(col = "none")
ggsave("../Manuscript/Figures/BIC_sigma_data.png", width = 4, height = 4)
