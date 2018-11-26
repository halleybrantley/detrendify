library(devtools)
load_all("detrendr")
rm(list=ls())
tau <- c(0.05, 0.1, 0.5, .9, 0.95)
nSim <- 100
n <- 300
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")

for(simDesign in simDesigns){                
  for (i in 1:nSim){
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    lam_SIC <- lambda_SIC(df$y, tau, 3)
    trend <- gurobi_trend(df$y, tau, lam_SIC$lambda, k=3)
    save(trend, lam_SIC,
         file = sprintf("../SimResults/detrend_SIC/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
    lam_valid <- lambda_valid(df$y, tau, 3, 5)
    trend <- gurobi_trend(df$y, tau, lam_valid$lambda, k=3)
    save(trend, lam_valid,
         file = sprintf("../SimResults/detrend_valid/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
  }
}
