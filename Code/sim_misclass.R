################################################################################
# Calculate mis-classification rates in peaks simulation
# Halley Brantley
################################################################################

rm(list=ls())
simDesign <- "peaks"
tau <- c(0.01, 0.05, 0.1)
nSim <- 100
methods <- c("detrend_eBIC", "detrend_SIC", "detrend_valid", "qsreg", "rqss", "npqw") 
miss_class <- as.data.frame(matrix(NA, nrow = nSim*length(methods), 
                             ncol = length(tau)+3))
thresh <- c(.8, 1, 1.2)
colnames(miss_class) <- c("Sim", "Method", "n", paste0("thresh_", thresh))
k <- 1

get_misclass <- function(trend, y, thresh, signal){
  y_adj <- y - trend
  signal_hat <- as.numeric(y_adj > thresh)
  return(mean(abs(signal-signal_hat)))
}


for (n in c(300,500,1000,5000)){
  for (i in 1:nSim){
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    df$signal <- as.numeric(df$peaks > 0.1)
    
    for (method in methods){
      load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
                   method, simDesign, n, i))

      miss_class[k,1] <- i
      miss_class[k,2] <- method
      miss_class[k,3] <- n
      miss_class[k,4] <- get_misclass(trend[,3], df$y, thresh[1], df$signal)
      miss_class[k,5] <- get_misclass(trend[,3], df$y, thresh[2], df$signal)
      miss_class[k,6] <- get_misclass(trend[,3], df$y, thresh[3], df$signal)
      k <- k+1
    }
  }
}


peaks_long <- miss_class %>% gather("tau", "missclass", -c("Sim", "Method", "n")) 

summary_peaks <- 
  peaks_long %>% group_by(Method, tau, n) %>% 
  summarise(
    mean_missclass = mean(missclass), 
    sd_missclass = sd(missclass)/sqrt(nSim)
  ) %>%
  ungroup() 


summary_peaks %>% 
  ggplot( aes(x = factor(n), y = mean_missclass, col = Method)) + 
  geom_point(position = position_dodge(width = 0.5)) +
  geom_linerange(aes(ymin = mean_missclass - 2*sd_missclass, ymax = mean_missclass + 2*sd_missclass), 
                 position = position_dodge(width = 0.5))+
  facet_grid(.~factor(tau))+
  theme_bw() +
  scale_color_brewer(palette = "Paired") +
  ylim(c(0.09, .19)) +
  labs(x = "", y="Miss-classified rate",  col = "Method") 

ggsave("../Manuscript/Figures/peaks_missclass.png", width = 10, height = 3)


n <- 500
i <- 75
method <- "detrend_eBIC"
load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
load(sprintf("../SimResults/%s/%s_n_%i_sim%03.0f.RData", 
             method, simDesign, n, i))
df$signal <- as.numeric(df$peaks > 0.1)
y_adj <- df$y - trend[,3]
signal_hat <- as.numeric(y_adj > thresh[2])

ggplot(df, aes(x=x, y=y)) + geom_line() +
  geom_point(data=subset(df, signal==1), col="red", size = 2) +
  geom_point(data=df[signal_hat==1,], col="blue")
ggsave("../Manuscript/Figures/peaks_eg_class.png", width = 10, height = 3)

