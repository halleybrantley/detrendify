library(tidyverse)
library(devtools)
library(jcolors)
library(Cairo)
library(grid)
rm(list=ls())
source("sim_generating_functions.R")
load_all("detrendr")
colPal <- c('#006d2c', '#2ca25f', '#66c2a4', 
            "#c2a5cf", "#9970ab", "#762a83")
text_size <- 14
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
nSim <- 100
simDesigns <- c( "mixednorm", "shapebeta", "gaus")
methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "rqss", "npqw", "qsreg")
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

p1 <- summary_stats %>% 
  filter(Design == "gaus") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_grid(.~factor(tau), scales = "free")+
  theme_bw() +
  theme(text = element_text(size=text_size), 
        axis.text.x = element_text(size=(text_size-5)),
        axis.title.y = element_text(margin = 
                                      margin(t = 0, r = 5, b = 0, l = 0)),
        plot.title = element_text(size = text_size))+
  scale_color_manual(values=colPal, breaks = methods) +
  labs(x = "n", y="RMSE", title = "Gaussian")+
  guides(col="none")+
  ylim(c(0,.153))

p2 <- summary_stats %>% 
  filter(Design == "shapebeta") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5)) +
  facet_grid(.~factor(tau), scales = "free")+
  theme_bw() +
  theme(text = element_text(size=text_size), 
        axis.text.x = element_text(size=(text_size-5)),
        plot.title = element_text(size = text_size)) +
  scale_color_manual(values=colPal, breaks = methods) +
  labs(x = "n", y="RMSE", title = "Beta") +
  guides(col="none")+
  ylim(c(0,.091))

p3 <- summary_stats %>% 
  filter(Design == "mixednorm") %>%
  ggplot( aes(x = factor(n), y = mean_mse, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_mse - 2*sd_mse, ymax = mean_mse + 2*sd_mse), 
                 position = position_dodge(width = 0.5))+
  facet_grid(.~factor(tau), scales = "free")+
  theme_bw() +
  theme(text = element_text(size=text_size), 
        axis.text.x = element_text(size=(text_size-5)),
        axis.title.y = element_text(margin = 
                                      margin(t = 0, r = 5, b = 0, l = 10)),
        plot.title = element_text(size = text_size), 
        legend.position = "bottom")+
  guides(col=guide_legend(nrow=2, byrow=TRUE)) + 
  scale_color_manual(values=colPal, breaks = methods) +
  labs(x = "n", y="RMSE", title = "Mixed Normal") + 
  ylim(c(0, .4))

fig.layout<-grid.layout(nrow=3,ncol=1, 
                        heights=c(2.8,2.8,4),  
                        widths=c(7,7,7), default.units="null", 
                        just=c("left","bottom"))

Cairo(file="../Manuscript/Figures/sim_metrics.png", 
      type="png",
      dpi = 400, 
      unit = "in",
      width=7, height=8)

pushViewport(viewport(layout=fig.layout))
print(p1, vp=viewport(layout.pos.row=1, layout.pos.col=1))
print(p2, vp=viewport(layout.pos.row=2,layout.pos.col=1))
print(p3, vp=viewport(layout.pos.row=3, layout.pos.col=1))
dev.off()
