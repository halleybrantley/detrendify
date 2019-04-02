library(fields)
rm(list=ls())
tau <- c(0.01, 0.05, 0.25, 0.5, .75, 0.95, 0.99)
tau <- c(0.01, 0.05, 0.1)
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")
simDesigns <- "peaks"
for (i in 1:100){
for (n in c(500,1000, 2000, 4000)){
for(simDesign in simDesigns){                
    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    trend <- matrix(NA, n, length(tau))
    for (j in 1:length(tau)){
      fit_qsreg <- qsreg(df$x, df$y, maxit.cv = 50, alpha=tau[j], hmin = -38)
      trend[,j] <- predict(fit_qsreg)    
    }
    save(trend, 
         file = sprintf("../SimResults/qsreg/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
}
}
}