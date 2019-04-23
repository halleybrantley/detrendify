library(devtools)
load_all("detrendr")
rm(list=ls())
# tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
tau <- c(0.01, 0.05, 0.1)
# simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")
simDesigns <- "peaks" 
i = 0
for (n in c(500,1000,2000,4000)){
  lambdaSeq = n^seq(0,1.4,length.out=15) 
  for (simDesign in simDesigns){                
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    trend <- matrix(NA, n, length(tau))
    trend_BIC <- list()
    for (j in 1:length(tau)){
      trend_BIC[[j]] <- get_trend_BIC(df$y, tau[j], 3, 
                            lambdaSeq = lambdaSeq, 
                            plot_lambda = FALSE, 
                            criteria = "eBIC")
      trend[,j] <- trend_BIC[[j]]$trend 
    }
    save(trend, trend_BIC,
         file = sprintf("../SimResults/detrend_cross/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
  }
}
