library(fields)
rm(list=ls())
tau <- c(0.05, 0.1, 0.5, .9, 0.95)
nSim <- 2
n <- 300
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")

for(simDesign in simDesigns){                
  for (i in 1:nSim){
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    trend <- matrix(NA, n, length(tau))
    for (j in 1:length(tau)){
      fit_qsreg <- qsreg(df$x, df$y, maxit.cv = 50, alpha=tau[j])
      trend[,j] <- predict(fit_qsreg)    
    }
    save(trend, 
         file = sprintf("../SimResults/qsreg/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
  }
}
