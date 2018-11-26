## $Id: example_npqw_numeric.R,v 1.4 2014/10/09 19:35:04 jracine Exp jracine $

source("lib_npqw_mixed.R")

regtype <- "ll"
ckertype <- "epanechnikov"
nmulti <- 2

set.seed(42)
n <- 1000
sigma <- 0.5
tau <- 0.10
x <- sort(runif(n))
dgp <- sin(2*pi*x)
dgp <- dgp/sd(dgp)

## You could replace y and X with data from your application. If you
## did so you might want to just use calls to
## plot(model,plot.errors.method="bootstrap") etc.

y <- dgp + rnorm(n,sd=sigma)
X <- data.frame(x)

## The method requires a pilot estimator of the conditional standard
## deviation, and this function provides it using the method of Fan &
## Yao (1998)

sd.pilot <- sd.fan.yao(x=X,
                       y=y,
                       regtype=regtype,
                       ckertype=ckertype,
                       nmulti=nmulti)

## We make a call to the function npqw.optim which takes the data,
## pilot estimator, and tau in (0,1) and computes cross-validated
## values of delta and the bandwidths.  Note that the li-racine kernel
## is used for simplicity as categorical bandwidths will lie in [0,1]

optim.out <- npqw.optim(x=X,
                        y=y,
                        tau=tau,
                        nmulti=nmulti,
                        sd.pilot=sd.pilot,
                        ckertype=ckertype,
                        regtype=regtype)

delta <- optim.out$delta
bws <- optim.out$bws

print(delta)
print(bws)

Q <- qnorm(delta,mean=y,sd=sd.pilot)

## Call to liracine kernel is necessary as these are used in the
## functions above

model <- npreg(txdat=X,
               tydat=Q,
               bws=bws,
               regtype=regtype,
               okertype="liracine",
               ukertype="liracine",               
               ckertype=ckertype)

plot(x,y,
     cex=0.25,
     sub=paste("tau = ",tau,", n = ",n,sep=""),
     col="grey")

lines(x,fitted(model),lwd=2)

