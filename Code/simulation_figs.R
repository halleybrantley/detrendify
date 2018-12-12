library(tidyverse)
library(devtools)
library(jcolors)
rm(list=ls())
source("trueQuantile.R")
load_all("detrendr")

tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
n <- 300
nSim <- 100
simDesigns <- c( "mixednorm", "shapebeta", "gaus")
methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "npqw", "qsreg", "rqss")
#methods <- "npqw"
#simDesigns <- "mixednorm"
MSEs <- as.data.frame(matrix(NA, nrow = nSim*length(methods)*length(simDesigns), 
                             ncol = length(tau)+4))
colnames(MSEs) <- c("Design", "Sim", "Method", "n", paste0("tau_", tau))
k <- 1
for(simDesign in simDesigns){  
  for (i in 1:nSim){
    for (n in c(300,500,1000)){
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    Q <- trueQuantile(simDesign, df$x, tau)
    for (method in methods){
      load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
              method, simDesign, n, i))
      MSEs[k,1] <- simDesign
      MSEs[k,2] <- i
      MSEs[k,3] <- method
      MSEs[k,4] <- n
      MSEs[k,5:ncol(MSEs)] <- colMeans((trend - Q)^2)
      k <- k+1
    }
  }
  }
}

tmp <- MSEs %>% filter(Design == "gaus", n==300)

hist(tmp[tmp$Method == "qsreg", "tau_0.5"])
median(tmp[tmp$Method == "npqw", "tau_0.5"])

which(tmp[tmp$Method == "npqw", "tau_0.5"] >.1)


MSEs_long <- MSEs %>% gather("tau", "MSE", -c("Design", "Sim", "Method", "n")) 
MSEs_long$RMSE <- sqrt(MSEs_long$MSE)

summary_stats <- 
  MSEs_long %>% group_by(Method, tau, Design, n) %>% 
  summarise(
    mean_mse = mean(RMSE), 
    sd_mse = sd(RMSE)/sqrt(nSim), 
    median_mse = median(RMSE), 
    mad_mse = median(abs(RMSE - median_mse))*1.482
  ) %>%
  ungroup() %>%
  mutate(tau_fac = tau, 
         tau = as.numeric(substr(tau_fac, 5, 10)))



summary_stats <- summary_stats %>% filter( tau > 0.01 & tau < 0.99) 

summary_stats %>% 
  filter(Design == "gaus") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_wrap(factor(tau)~., scales = "free", ncol=5)+
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  labs(x = "n", y="RMSE", title = "Gaussian")
ggsave("../Manuscript/Figures/gaus_mse.png", width = 10, height = 3)

summary_stats %>% 
  filter(Design == "shapebeta") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5)) +
  facet_wrap(factor(tau)~., scales = "free", ncol=5)+
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  labs(x = "n", y="RMSE", title = "Beta")
ggsave("../Manuscript/Figures/shapebeta_mse.png", width = 10, height = 3)

summary_stats %>% 
  filter(Design == "mixednorm") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_wrap(factor(tau)~., scales = "free", ncol=5)+
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  labs(x = "n", y="RMSE", title = "Mixed Normal")
ggsave("../Manuscript/Figures/mixednorm_mse.png", width = 10, height = 3)


################################################################################
i <- 10
n <- 500
tau <- c(0.01, 0.05, 0.1)
simDesign <- "peaks"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             "detrend_eBIC", simDesign, n, i))
trend_detrend <- as.data.frame(trend) 
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             "qsreg", simDesign, n, i))
trend_qsreg <- as.data.frame(trend) 

trend_detrend$x <- df$x
trend_detrend$q <- df$baseline+qnorm(tau[3], sd = 0.25)
trend_detrend$q2 <- df$baseline+qnorm(tau[2], sd = 0.25)
trend_detrend$qsreg <- trend_qsreg[,3]
trend_detrend$qsreg2 <- trend_qsreg[,2]

ggplot(df, aes(x=x, y=y)) + 
  geom_line(col="grey") +
  geom_line(data=trend_detrend, aes(y=q, x=x, col = "true 0.1"), lwd = 1.5) +
  geom_line(data=trend_detrend, aes(y=q2, x=x, col = "true 0.05"), lwd = 1.5) +
  geom_line(data=trend_detrend, aes(y=V3, x=x, col = "detrendr 0.1")) +
  geom_line(data=trend_detrend, aes(y=qsreg, x=x, col = "qsreg 0.1")) + 
  geom_line(data=trend_detrend, aes(y=V2, x=x, col = "detrendr 0.05")) +
  geom_line(data=trend_detrend, aes(y=qsreg2, x=x, col = "qsreg 0.05")) +
  theme_bw() + 
  scale_color_brewer(palette = "Set1") +
  labs(col="")
ggsave("../Manuscript/Figures/ex_baseline.png", width = 7, height = 4)

mean((trend_detrend$V3 - trend_detrend$q)^2)  
mean((trend_detrend$qsreg - trend_detrend$q)^2)  

################################################################################
simDesign <- "peaks"
tau <- c(0.01, 0.05, 0.1)
nSim <- 100
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
plot(y~x, df, col="grey", type="l")
lines(baseline~x, df, col="red")
lines((peaks+baseline)~x, df, col="blue")

methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "qsreg", "rqss", "npqw") 
MSEs <- as.data.frame(matrix(NA, nrow = nSim*length(methods), 
                             ncol = length(tau)+3))
colnames(MSEs) <- c("Sim", "Method", "n", paste0("tau_", tau))
k <- 1
  for (n in c(300,500,1000,5000,10000)){
    for (i in 1:nSim){
      if (i == 49){ next }
      load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))

      for (method in methods){
        load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
                     method, simDesign, n, i))
        MSEs[k,1] <- i
        MSEs[k,2] <- method
        MSEs[k,3] <- n
        MSEs[k,4] <- mean((trend[,1] - (df$baseline+qnorm(tau[1], sd = 0.25)))^2)
        MSEs[k,5] <- mean((trend[,2] - (df$baseline+qnorm(tau[2], sd = 0.25)))^2)
        MSEs[k,6] <- mean((trend[,3] - (df$baseline+qnorm(tau[3], sd = 0.25)))^2)
        k <- k+1
      }
    }
  }

tmp <- MSEs %>% filter(n<10000)

hist(tmp[tmp$Method == "detrend_eBIC", "tau_0.05"], 50)
hist(tmp[tmp$Method == "detrend_eBIC", "tau_0.05"], 50, col="blue", add=T)

plot(tmp[tmp$Method == "detrend_eBIC", "tau_0.1"]~
       tmp[tmp$Method == "detrend_SIC", "tau_0.1"])
abline(0,1)
which(tmp[tmp$Method == "detrend_SIC", "tau_0.05"] > .6)

peaks_long <- MSEs %>% gather("tau", "MSE", -c("Sim", "Method", "n")) %>%
  filter(n < 10000)
peaks_long$RMSE <- sqrt(peaks_long$MSE)
summary_peaks <- 
  peaks_long %>% group_by(Method, tau, n) %>% 
  summarise(
    mean_mse = mean(RMSE), 
    sd_mse = sd(RMSE)/sqrt(nSim), 
    median_mse = median(RMSE), 
    mad_mse = median(abs(RMSE - median_mse))*1.482
  ) %>%
  ungroup() %>%
  mutate(tau_fac = tau, 
         tau = as.numeric(substr(tau_fac, 5, 10)))

summary_peaks %>% 
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_wrap(factor(tau)~., scales = "free", ncol=5)+
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  labs(x = "", y="RMSE",  col = "Method")
ggsave("../Manuscript/Figures/peaks_mse.png", width = 10, height = 3)
