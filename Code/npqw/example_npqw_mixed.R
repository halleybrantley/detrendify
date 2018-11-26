## $Id: example_npqw_mixed.R,v 1.4 2014/10/09 19:34:57 jracine Exp jracine $

source("lib_npqw_mixed.R")

regtype <- "ll"
ckertype <- "gaussian"
nmulti <- 2

set.seed(42)
n <- 500
sigma <- 0.25
tau.seq <- c(0.01,0.25,0.5,0.75,0.99)
x1 <- runif(n)
x2 <- runif(n)
x3 <- factor(rbinom(n,1,.5))
x4 <- rnorm(n)
dgp <- sin(2*pi*x1) + 2*x2
dgp <- dgp/sd(dgp)

## You could replace y and X with data from your application. If you
## did so you might want to just use calls to
## plot(model,plot.errors.method="bootstrap") etc.

y <- dgp + rnorm(n,sd=sigma)
X <- data.frame(x1,x2,x3,x4)

## The method requires a pilot estimator of the conditional standard
## deviation, and this function provides it using the method of Fan &
## Yao (1998)

sd.pilot <- sd.fan.yao(x=X,
                       y=y,
                       regtype=regtype,
                       ckertype=ckertype,
                       nmulti=nmulti)

plot(x1,y,
     cex=0.25,
     col="grey")

for(i in 1:length(tau.seq)) {

    ## For each value of tau in tau.seq, we cross-validate bandwidths
    ## and delta

    optim.out <- npqw.optim(x=X,
                            y=y,
                            tau=tau.seq[i],
                            nmulti=nmulti,
                            sd.pilot=sd.pilot,
                            ckertype=ckertype,
                            regtype=regtype)


    ## Note that the liracine kernel is used for simplicity as
    ## categorical bandwidths will lie in [0,1]
    
    delta <- optim.out$delta
    bws <- optim.out$bws

    print(delta)
    print(bws)

    ## Here we construct the evaluation data for one predictor holding
    ## the remaining constant at their medians

    x1.eval <- seq(0,1,length=100)
    x2.eval <- median(x2)
    x3.eval <- factor(uocquantile(x3,0.5),levels=levels(x3))
    x4.eval <- median(x4)
    x.eval <- data.frame(x1=x1.eval,x2=x2.eval,x3=x3.eval,x4=x4.eval)

    Q <- qnorm(delta,mean=y,sd=sd.pilot)

    model <- npreg(txdat=X,
                   tydat=y,
                   exdat=x.eval,
                   bws=bws,
                   regtype=regtype,
                   ckertype=ckertype,
                   okertype="liracine",
                   ukertype="liracine")
    
    Q.dgp <- qnorm(tau.seq[i], mean=sin(2*pi*x1.eval) + 2*median(x2),sd=sigma)

    lines(x1.eval,fitted(model),col=i,lwd=2)
    lines(x1.eval,Q.dgp,lty=2,col="lightgrey",lwd=2)

}

legend("topleft",
       paste("tau = ",tau.seq,sep=""),
       lty=1,
       col=1:length(tau.seq),
       lwd=rep(2,length(tau.seq)),
       bty="n")

