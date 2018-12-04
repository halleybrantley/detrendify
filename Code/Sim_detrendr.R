library(devtools)
load_all("detrendr")
rm(list=ls())
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
nSim <- 100
i = 0
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")

for(simDesign in simDesigns){                
  for (n in c(300,500,1000)){
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    lam_SIC <- lambda_SIC(df$y, tau, 3, 
                          lambdaSeq = n^seq(1,2,length.out=10), 
                          plot_lambda = TRUE, 
                          single_lambda = TRUE)
    trend <- gurobi_trend(df$y, tau, lam_SIC$lambda, k=3)
    save(trend, lam_SIC,
         file = sprintf("../SimResults/detrend_SIC/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
    lam_valid <- lambda_valid(df$y, tau, 3, 5, lambdaSeq = n^seq(1,2,length.out=10),
                              single_lambda = TRUE)
    trend <- gurobi_trend(df$y, tau, lam_valid$lambda, k=3)
    save(trend, lam_valid,
         file = sprintf("../SimResults/detrend_valid/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
  }
}
