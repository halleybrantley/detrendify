library(devtools)
load_all("detrendr")
rm(list=ls())
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
tau <- c(0.01, 0.05, 0.1)
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")
simDesigns <- "peaks" 

i = 1
for (n in c(500,1000,2000,4000)){
  lambdaSeq = n^seq(0,1.4,length.out=15) 
  overlap <- n/10
  window_size <- round(n/3+overlap*(3-1)/3)
  for (simDesign in simDesigns){                
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    trend_BIC <- get_windows_BIC(df$y, tau, k=3,
                          lambdaSeq = lambdaSeq,
                          window_size = window_size, 
                          overlap = overlap,
                          plot_lambda = FALSE,
                          criteria = "eBIC")
    trend <- trend_BIC$trend
    save(trend, trend_BIC,
         file = sprintf("../SimResults/windows/%s_n_%i_sim%03.0f.RData",
                        simDesign, n, i))
  }
}
