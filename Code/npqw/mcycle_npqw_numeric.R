## $Id: mcycle_npqw_numeric.R,v 1.1 2014/10/10 15:40:18 jracine Exp jracine $

## Illustration for the mcycle data 

source("lib_npqw_mixed.R")
data(mcycle,package="MASS")
attach(mcycle)

set.seed(42)

nmulti <- 5
ckertype <- "epanechnikov"
regtype <- "ll"
tau.vec <- c(0.01,0.05,0.25,0.50,0.75,0.95,0.99)

require(np)
n <- nrow(mcycle)

plot(times,accel,
     main="",
     xlab="times",
     ylab="accel",
     ylim=c(-175,100),
     sub=paste("kernel = ",ckertype,", regtype = ",regtype,", nmulti = ",nmulti,sep=""))

bws.vec <- numeric(length(tau.vec))
delta.vec <- numeric(length(tau.vec))

sd.pilot <- sd.fan.yao(x=times,
                       y=accel,
                       nmulti=nmulti,
                       regtype=regtype,
                       ckertype=ckertype)
        
for(i in 1:length(tau.vec)) {

    optim.out <- npqw.optim(x=times,
                            y=accel,
                            tau=tau.vec[i],
                            nmulti=nmulti,
                            sd.pilot=sd.pilot,
                            ckertype=ckertype,
                            regtype=regtype)

    bws.vec[i] <- optim.out$bws
    delta.vec[i] <- optim.out$delta

    Q <- qnorm(delta.vec[i],mean=accel,sd=sd.pilot)

    print(paste("tau = ",tau.vec[i],
                ", best fval = ",formatC(optim.out$value,format="f",digits=4),
                ", bws = ",formatC(bws.vec[i],format="f",digits=4),
                ", delta = ",formatC(delta.vec[i],format="f",digits=4),sep=""))
    
    lines(times,fitted(npreg(Q~times,
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

