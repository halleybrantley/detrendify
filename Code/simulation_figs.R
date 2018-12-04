library(tidyverse)
rm(list=ls())
source("trueQuantile.R")

tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
n <- 300
nSim <- 100
simDesigns <- c("gaus", "mixednorm", "shapebeta")
methods <- c("detrend_SIC", "detrend_valid", "npqw", "qsreg", "rqss") #, "rqss")
MSEs <- as.data.frame(matrix(NA, nrow = nSim*length(methods)*length(simDesigns), 
                             ncol = length(tau)+4))
colnames(MSEs) <- c("Design", "Sim", "Method", "n", paste0("tau_", tau))
k <- 1
for(simDesign in simDesigns){  
  for (n in c(300,500,1000)){
  for (i in 1:nSim){
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

tmp <- MSEs %>% filter(Design == "gaus", n==1000)

hist(tmp[tmp$Method == "detrend_SIC", "tau_0.01"])
mean(tmp[tmp$Method == "rqss", "tau_0.01"])

which(tmp[tmp$Method == "detrend_SIC", "tau_0.01"]>.04)


MSEs_long <- MSEs %>% gather("tau", "MSE", -c("Design", "Sim", "Method", "n")) 
summary_stats <- 
  MSEs_long %>% group_by(Method, tau, Design, n) %>% 
  summarise(
    mean_mse = mean(MSE), 
    sd_mse = sd(MSE)/sqrt(nSim), 
    median_mse = median(MSE), 
    mad_mse = median(abs(MSE - median_mse))*1.482
  ) %>%
  ungroup() %>%
  mutate(tau_fac = tau, 
         tau = as.numeric(substr(tau_fac, 5, 10)))



ggplot(summary_stats, aes(x=n, y = median_mse*100, col = Method)) + 
  geom_point(position = position_dodge(width = 0.2)) +
  geom_linerange(aes(ymin = median_mse*100 - mad_mse*100, 
                     ymax = median_mse*100 + mad_mse*100), 
                 position = position_dodge(width = 0.2))+
  facet_grid(Design~tau_fac, scales = "free")+
  theme_bw() +
  scale_color_brewer(palette = "Set1") 

summary_stats %>% 
  filter(Design == "gaus") %>%
  ggplot( aes(x = tau_fac, y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_grid(n~., scales = "free")+
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "", y="MSE", title = "Gaussian")
ggsave("../Manuscript/Figures/gaus_mse.png", width = 6, height = 5)

summary_stats %>% 
  filter(Design == "shapebeta") %>%
  ggplot( aes(x = tau_fac, y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_grid(n~., scales = "free")+
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "", y="MSE", title = "Beta")
ggsave("../Manuscript/Figures/shapebeta_mse.png", width = 6, height = 5)

summary_stats %>% 
  filter(Design == "mixednorm") %>%
  ggplot( aes(x = tau_fac, y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_grid(n~., scales = "free")+
  theme_bw() +
  scale_color_brewer(palette = "Set1") +
  labs(x = "", y="MSE", title = "Mixed Normal")
ggsave("../Manuscript/Figures/mixednorm_mse.png", width = 6, height = 5)


i <- 1
simDesign <- "gaus"
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             "detrend_SIC", simDesign, n, i))
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
Q <- trueQuantile(simDesign, df$x, tau)
trend_detrend <- as.data.frame(trend) 
trend_detrend$x <- df$x
ggplot(df, aes(x=x, y=y)) + 
  geom_line(col="grey") +
  geom_line(data=trend_rqss, aes(y=V1, x=x))+
  geom_line(data=trend_detrend, aes(y=V1, x=x))
  


load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
plot(y~x, df, col="grey", type="l")
lines(baseline~x, df, col="red")
lines((peaks+baseline)~x, df, col="blue")

method <- "npqw"
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             method, simDesign, n, i))
lines(trend[,1]~df$x)
