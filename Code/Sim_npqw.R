
rm(list=ls())
source("npqw/lib_npqw_mixed.R")

tau <- c(0.05, 0.1, 0.5, .9, 0.95)
nSim <- 2
n <- 300
simDesigns <- c("gaus", "mixednorm", "shapebeta", "peaks")

regtype <- "ll"
ckertype <- "epanechnikov"
nmulti <- 2

i = 0 
################################################################################
  for(simDesign in simDesigns){                

    load(sprintf("../SimData/%s_n_%i_sim%03.0f.RData", simDesign, n, i))
    trend <- matrix(NA, n, length(tau))
    sd.pilot <- sd.fan.yao(x=df$x,
                           y=df$y,
                           regtype=regtype,
                           ckertype=ckertype,
                           nmulti=nmulti)
    
    ## We make a call to the function npqw.optim which takes the data,
    ## pilot estimator, and tau in (0,1) and computes cross-validated
    ## values of delta and the bandwidths.  Note that the li-racine kernel
    ## is used for simplicity as categorical bandwidths will lie in [0,1]
    

    
    for (j in 1:length(tau)){

      optim.out <- npqw.optim(x=df$x,
                              y=df$y,
                              tau=tau[j],
                              nmulti=nmulti,
                              sd.pilot=sd.pilot,
                              ckertype=ckertype,
                              regtype=regtype)
      delta <- optim.out$delta
      bws <- optim.out$bws
      Q <- qnorm(delta,mean=df$y,sd=sd.pilot)
      
      ## Call to liracine kernel is necessary as these are used in the
      ## functions above
      
      model <- npreg(txdat=df$x,
                     tydat=Q,
                     bws=bws,
                     regtype=regtype,
                     okertype="liracine",
                     ukertype="liracine",               
                     ckertype=ckertype)
      
      trend[,j] <- fitted(model)
    }
    
    save(trend, 
         file = sprintf("../SimResults/npqw/%s_n_%i_sim%03.0f.RData", 
                        simDesign, n, i))
  }
################################################################################
