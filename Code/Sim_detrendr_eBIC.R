library(devtools)
load_all("detrendr")
rm(list=ls())
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
# nSim <- 100
# n <- 500
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")
i = 0
for (n in c(300, 500, 1000)){
for (simDesign in simDesigns){
    if (simDesign == "peaks"){
      tau <- c(0.01, 0.05, 0.1)
    }
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    lam <- lambda_eBIC(df$y, tau, 3,  gamma = 0.1,
                          plot_lambda = FALSE, single_lambda = FALSE)
    trend <- gurobi_trend(df$y, tau, lam$lambda, k=3)
    save(trend, lam,
         file = sprintf("../SimResults/detrend_eBIC/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
  }
}
