library(quantreg)
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
    trend <- matrix(NA, n, length(tau))
    for (j in 1:length(tau)){
      lam_SIC <- lambda_SIC(df$y, tau[j], 2, 
                             lambdaSeq = seq(5, n/2, n/25), 
                            plot_lambda = FALSE)
      fit_rqss <- rqss(y ~ qss(x, lambda = 2*lam_SIC$lambda/n), tau = tau[j], data = df)
      trend[,j] <- predict(fit_rqss, df)  
    }
    save(trend, 
         file = sprintf("../SimResults/rqss/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
  }
}
