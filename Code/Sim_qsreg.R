library(fields)
rm(list=ls())
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
nSim <- 2
n <- 600
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")
i = 0

for(simDesign in simDesigns){                
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
