library(tidyverse)
library(devtools)
library(jcolors)
rm(list=ls())
source("sim_generating_functions.R")
load_all("detrendr")
colPal <- rev(c('#762a83','#9970ab','#c2a5cf',
                '#a6dba0','#5aae61','#1b7837'))

tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
nSim <- 100
simDesigns <- c( "mixednorm", "shapebeta", "gaus")
methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "qsreg", "rqss", "npqw")
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
         tau = as.numeric(substr(tau_fac, 5, 10)), 
         Method = factor(Method, levels = methods))



summary_stats <- summary_stats %>% filter( tau > 0.01 & tau < 0.99) 

summary_stats %>% 
  filter(Design == "gaus") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_grid(.~factor(tau), scales = "free")+
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  labs(x = "n", y="RMSE", title = "Gaussian")+
  ylim(c(0,.153))
ggsave("../Manuscript/Figures/gaus_mse.png", width = 10, height = 3)

summary_stats %>% 
  filter(Design == "shapebeta") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5)) +
  facet_grid(.~factor(tau), scales = "free")+
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  labs(x = "n", y="RMSE", title = "Beta") +
  ylim(c(0,.091))
ggsave("../Manuscript/Figures/shapebeta_mse.png", width = 10, height = 3)

summary_stats %>% 
  filter(Design == "mixednorm") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_grid(.~factor(tau), scales = "free")+
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  labs(x = "n", y="RMSE", title = "Mixed Normal") + 
  ylim(c(0, .4))
ggsave("../Manuscript/Figures/mixednorm_mse.png", width = 10, height = 3)


