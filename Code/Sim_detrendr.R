library(devtools)
load_all("detrendr")
rm(list=ls())
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
tau <- c(0.01, 0.05, 0.1)
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")
simDesigns <- "peaks" 
i = 0
for (n in c(500,1000,2000,4000)){
  lambdaSeq = n^seq(0,1.4,length.out=15) 
  for (simDesign in simDesigns){                
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    # trend_BIC <- get_trend_BIC(df$y, tau, 3, 
    #                       lambdaSeq = lambdaSeq, 
    #                       plot_lambda = FALSE, 
    #                       criteria = "eBIC")
    # trend <- trend_BIC$theta 
    # save(trend, trend_BIC,
    #      file = sprintf("../SimResults/detrend_eBIC/%s_n_%i_sim%03.0f.RData", 
    #                     simDesign, n, i))
    trend_valid <- get_trend_BIC(y=df$y, tau, k=3, 
                                 lambdaSeq = lambdaSeq,
                                 plot_lambda = TRUE, 
                                 criteria = "valid")
    trend <- trend_valid$theta
    save(trend, trend_valid,
         file = sprintf("../SimResults/detrend_valid/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
    
    # trend_SIC <- get_trend_BIC(df$y, tau, 3, 
    #                            lambdaSeq = lambdaSeq,
    #                            plot_lambda = FALSE, 
    #                            criteria = "SIC")
    # trend <- trend_SIC$theta
    # save(trend, trend_SIC,
    #      file = sprintf("../SimResults/detrend_SIC/%s_n_%i_sim%03.0f.RData", 
    #                     simDesign, n, i))
  }
}
