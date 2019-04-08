################################################################################
# Peaks Simulation Results
# Halley Brantley
################################################################################
library(tidyverse)
rm(list=ls())
i <- 10
n <- 4000
tau <- c(0.01, 0.05, 0.1)
simDesign <- "peaks"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             "detrend_valid", simDesign, n, i))
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
  #geom_line(data=trend_detrend, aes(y=q, x=x, col = "true 0.1"), lwd = 1.5) +
  geom_line(data=trend_detrend, aes(y=q2, x=x, col = "true 0.05"), lwd = 1.5) +
  #geom_line(data=trend_detrend, aes(y=V3, x=x, col = "detrendr 0.1")) +
  #geom_line(data=trend_detrend, aes(y=qsreg, x=x, col = "qsreg 0.1")) + 
  geom_line(data=trend_detrend, aes(y=V2, x=x, col = "detrendr 0.05")) +
  geom_line(data=trend_detrend, aes(y=qsreg2, x=x, col = "qsreg 0.05")) +
  theme_bw() + 
  scale_color_brewer(palette = "Set1") +
  labs(col="")
ggsave("../Manuscript/Figures/ex_baseline.png", width = 7, height = 4)

################################################################################
simDesign <- "peaks"
tau <- c(0.01, 0.05, 0.1)
nSim <- 100
colPal <- rev(c('#762a83','#9970ab','#c2a5cf', '#d9f0d3',
                '#a6dba0','#5aae61','#1b7837'))
methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "detrend_Xing",
             "rqss", "npqw", "qsreg") 
MSEs <- as.data.frame(matrix(NA, nrow = nSim*length(methods), 
                             ncol = length(tau)+3))
colnames(MSEs) <- c("Sim", "Method", "n", paste0("tau_", tau))
k <- 1
for (n in c(500,1000,2000,4000)){
  for (i in 1:nSim){
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


peaks_long <- MSEs %>% gather("tau", "MSE", -c("Sim", "Method", "n")) 
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
         tau = as.numeric(substr(tau_fac, 5, 10)), 
         Method = factor(Method, levels = methods))

text_size <- 14

summary_peaks %>% 
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_grid(.~factor(tau), scales = "free")+
  theme_bw() +
  scale_color_manual(values=colPal, breaks = methods) +
  theme(text = element_text(size=text_size), 
        axis.text.x = element_text(size=(text_size-5)),
        plot.title = element_text(size = text_size)) +
  labs(x = "", y="RMSE",  col = "Method") +
  ylim(c(0,1.2))
ggsave("../Manuscript/Figures/peaks_mse.png", width = 7, height = 2.5)



