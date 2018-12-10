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
    lam <- lambda_eBIC(df$y, tau, 3,  gamma = 1,
                       lambdaSeq = n^seq(0, 1.1, length.out=50),
                       plot_lambda = TRUE, single_lambda = FALSE)
    trend <- gurobi_trend(df$y, tau, lam$lambda, k=3)
    save(trend, lam,
         file = sprintf("../SimResults/detrend_eBIC/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
  }
}

plot(y~x, df, type="l", col="grey")
for(i in 1:length(tau)){
  lines(trend[,i]~df$x)
}
