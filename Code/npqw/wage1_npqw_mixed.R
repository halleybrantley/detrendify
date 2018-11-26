## $Id: wage1_npqw_mixed.R,v 1.2 2014/10/09 19:33:05 jracine Exp jracine $

## $Id: wage1_npqw_mixed.R,v 1.2 2014/10/09 19:33:05 jracine Exp jracine $

source("lib_npqw_mixed.R")

pdf(file="wage1.pdf")

regtype <- "ll"
ckertype <- "gaussian"
nmulti <- 5

tau.seq <- c(0.05,0.25,0.5,0.75,0.95)

data(wage1)
attach(wage1)

y <- lwage
X <- data.frame(educ,exper,tenure,female)

## The method requires a pilot estimator of the conditional standard
## deviation, and this function provides it using the method of Fan &
## Yao (1998)

sd.pilot <- sd.fan.yao(x=X,
                       y=y,
                       regtype=regtype,
                       ckertype=ckertype,
                       nmulti=nmulti)

for(i in 1:length(tau.seq)) {

    ## We make a call to the function npqw.optim which takes the data,
    ## pilot estimator, and tau in (0,1) and computes cross-validated
    ## values of delta and the bandwidths.  Note that the liracine
    ## kernel is used for simplicity as categorical bandwidths will
    ## lie in [0,1]

    optim.out <- npqw.optim(x=X,
                            y=y,
                            tau=tau.seq[i],
                            nmulti=nmulti,
                            sd.pilot=sd.pilot,
                            ckertype=ckertype,
                            regtype=regtype)
    
    bws <- optim.out$bws
    delta <- optim.out$delta

    print(delta)
    print(bws)

    Q <- qnorm(delta,mean=y,sd=sd.pilot)

    ## Call to liracine kernel in npreg is necessary as these are used
    ## in the functions above

    model <- npreg(txdat=X,
                   tydat=Q,
                   bws=bws,
                   regtype=regtype,
                   okertype="liracine",
                   ukertype="liracine",                        
                   ckertype=ckertype)

    summary(model)
    summary(model$bws)    

    plot(model,plot.errors.method="bootstrap",plot.errors.boot.num=50)
    plot(model,gradients=TRUE,plot.errors.method="bootstrap",plot.errors.boot.num=50)

}
