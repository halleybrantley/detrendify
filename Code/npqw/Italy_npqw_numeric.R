## $Id: Italy_npqw_numeric.R,v 1.2 2014/10/10 17:28:12 jracine Exp jracine $

## Illustration for the Italian GDP data

source("lib_npqw_mixed.R")
require(np)
data(Italy)
n <- nrow(Italy)
attach(Italy)

nmulti <- 5
ckertype <- "epanechnikov"
regtype <- "ll"
tau.vec <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)

plot(year,
     gdp,
     main="",
     xlab="year",
     ylab="gdp",
     sub=paste("kernel = ",ckertype,", regtype = ",regtype,", nmulti = ",nmulti,sep=""))

year <- as.numeric(as.character(year))

bws.vec <- numeric(length(tau.vec))
delta.vec <- numeric(length(tau.vec))

sd.pilot <- sd.fan.yao(x=year,
                       y=gdp,
                       nmulti=nmulti,
                       regtype=regtype,
                       ckertype=ckertype)
        
for(i in 1:length(tau.vec)) {

    optim.out <- npqw.optim(x=year,
                            y=gdp,
                            tau=tau.vec[i],
                            nmulti=nmulti,
                            sd.pilot=sd.pilot,
                            ckertype=ckertype,
                            regtype=regtype)

    bws.vec[i] <- optim.out$bws
    delta.vec[i] <- optim.out$delta

    Q <- qnorm(delta.vec[i],mean=gdp,sd=sd.pilot)

    print(paste("tau = ",tau.vec[i],
                ", best fval = ",formatC(optim.out$value,format="f",digits=4),
                ", bws = ",formatC(bws.vec[i],format="f",digits=4),
                ", delta = ",formatC(delta.vec[i],format="f",digits=4),sep=""))
    
    lines(ordered(year),fitted(npreg(Q~year,
                                     bws=bws.vec[i],
                                     ckertype=ckertype,
                                     okertype="liracine",
                                     ukertype="liracine",
                                     regtype=regtype)),
          col=i,
          lty=i)

}

legend("topleft",
       paste("tau = ",tau.vec,
             " (h = ",formatC(bws.vec,format="f",digits=3),
             ", delta = ",formatC(delta.vec,format="f",digits=3),")",sep=""),
       col=1:i,
       lty=1:i,
       bty="n",
       cex=0.75)

