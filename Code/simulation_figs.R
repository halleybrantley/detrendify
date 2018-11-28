library(tidyverse)
rm(list=ls())
source("trueQuantile.R")

tau <- c(0.05, 0.1, 0.5, .9, 0.95)
n <- 300
nSim <- 100
simDesigns <- c("gaus", "mixednorm", "shapebeta")
x <- seq(0.5, n, 1)/n
methods <- c("detrend_SIC", "detrend_valid", "npqw", "qsreg", "rqss") #, "rqss")
MSEs <- as.data.frame(matrix(NA, nrow = nSim*length(methods)*length(simDesigns), 
                             ncol = length(tau)+3))
colnames(MSEs) <- c("Design", "Sim", "Method", paste0("tau_", tau))
k <- 1
for(simDesign in simDesigns){  
  Q <- trueQuantile(simDesign, x, tau)
  for (i in 1:nSim){
    for (method in methods){
      load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
              method, simDesign, n, i))
      MSEs[k,1] <- simDesign
      MSEs[k,2] <- i
      MSEs[k,3] <- method
      MSEs[k,4:ncol(MSEs)] <- colMeans((trend - Q)^2)
      k <- k+1
    }
  }
}

MSEs_long <- MSEs %>% gather("tau", "MSE", -c("Design", "Sim", "Method")) 
summary_stats <- 
  MSEs_long %>% group_by(Method, tau, Design) %>% 
  summarise(
    mean_mse = mean(MSE), 
    sd_mse = sd(MSE)/sqrt(nSim), 
    median_mse = median(MSE), 
    mad_mse = median(abs(MSE - median_mse))*1.482
  ) %>%
  ungroup() %>%
  mutate(tau_fac = tau, 
         tau = as.numeric(substr(tau_fac, 5, 10)))

tmp <- MSEs_long %>%
  filter(tau == "tau_0.05", Method == "npqw", Design == "gaus") %>%
  select(MSE) 
  
hist(tmp$MSE, 20)
median(tmp$MSE[-62])*100

summary_stats %>%
  filter(tau_fac == "tau_0.95", Method == "npqw", Design == "gaus") %>%
  select(median_mse, mad_mse)*100


ggplot(summary_stats, aes(x =tau_fac, y = median_mse*100, col = Method)) + 
  geom_point(position = position_dodge(width = 0.2)) +
  geom_linerange(aes(ymin = median_mse*100 - mad_mse*100, 
                     ymax = median_mse*100 + mad_mse*100), 
                 position = position_dodge(width = 0.2))+
  facet_grid(Design~., scales = "free")+
  theme_bw() +
  scale_color_brewer(palette = "Set1") 

ggplot(summary_stats, aes(x = tau_fac, y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.2)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.2))+
  facet_grid(Design~., scales = "free")+
  theme_bw() +
  scale_color_brewer(palette = "Set1")
