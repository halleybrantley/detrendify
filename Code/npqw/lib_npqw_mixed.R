## $Id: lib_npqw_mixed.R,v 1.4 2014/10/09 20:36:32 jracine Exp jracine $

## Proof of concept code for multivariate direct quantile weighted
## kernel estimation with mixed covariate types.

require(np)
options(np.tree=TRUE)

## Fan-Yao (1998) estimator of conditional standard deviation of Y -
## fit conditional mean, regress squared residuals on
## predictor. Fourth order kernel used for asymptotic reasons.

## NOTE: 08/04/14 - np version 0.60-2 appears to suffer from a glitch
## in local linear regression with higher order kernels, Tristen has
## been alerted avoid use of higher order kernels with the pilot
## estimator for the moment being.

sd.fan.yao <- function(x,
                       y,
                       nmulti=1,
                       regtype="lc",
                       ckertype="gaussian",
                       ckerorder=2) {

    if(missing(x)) stop("must provide x")
    if(missing(y)) stop("must provide y")

    r.xy <- residuals(npreg(tydat=y,
                            txdat=x,
                            ckertype=ckertype,
                            regtype=regtype,
                            ckerorder=ckerorder,
                            okertype="liracine",
                            ukertype="liracine",
                            nmulti=nmulti))

    fan.yao <- fitted(npreg(tydat=r.xy**2,
                            txdat=x,
                            ckertype=ckertype,
                            regtype=regtype,
                            ckerorder=ckerorder,
                            okertype="liracine",
                            ukertype="liracine",
                            nmulti=nmulti))

    return(sqrt(abs(fan.yao)))

}

NZD <- function(a) {
    sapply(1:NROW(a), function(i) {if(a[i] < 0) min(-.Machine$double.eps,a[i]) else max(.Machine$double.eps,a[i])})
}

## Check function

check.function <- function(u,tau=0.5) {

#    if(missing(u)) stop(" Error: u must be provided")
#    if(tau <= 0 | tau >= 1) stop(" Error: tau must lie in (0,1)")

    return(u*(tau-ifelse(u<0,1,0)))

}

##  Delete-one direct weighted quantile estimator

npqw.estimator.loo <- function(params,
                               x,
                               y,
                               tau=0.5,
                               sd.pilot,
                               ckertype="gaussian",
                               regtype="lc") {

    if(missing(x)) stop("must provide x")
    if(missing(y)) stop("must provide y")
    if(missing(params)) stop("must provide params")
    if(missing(sd.pilot)) stop("must provide sd.pilot")

    delta <- params[1]
    bws <- params[-1]

    Q <- qnorm(delta,mean=y,sd=sd.pilot)

    cv.maxPenalty <- sqrt(.Machine$double.xmax)

    x.col.numeric <- sapply(1:ncol(x),function(i){is.numeric(x[,i])})
    k <- ncol(as.data.frame(x[,x.col.numeric]))

    if(regtype=="lc") {
        p <- 0
    } else {
        p <- 1
    }

    ## With no numeric predictors local linear is local constant

    if(regtype=="lc" || k==0) {

        ## Local constant estimation 

        tww <- npksum(txdat = x,
                      weights = as.matrix(data.frame(1,Q)),
                      tydat = rep(1,length(Q)),
                      bws = bws,
                      bandwidth.divide = TRUE,
                      ckertype=ckertype,
                      okertype="liracine",
                      ukertype="liracine",
                      leave.one.out = TRUE)$ksum

        m <- tww[2,]/NZD(tww[1,])

    } else {

        ## Local linear estimation 

        n <- length(y)

        W <- as.matrix(cbind(1,x[,x.col.numeric]))

        tww <- npksum(txdat = x,
                      tydat = as.matrix(cbind(Q,W)),
                      weights = W,
                      bws = bws,
                      bandwidth.divide = TRUE,
                      ckertype=ckertype,
                      okertype="liracine",
                      ukertype="liracine",
                      leave.one.out = TRUE)$ksum

        tyw <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[1,,]
        tww <- array(tww,dim = c(ncol(W)+1,ncol(W),n))[-1,,]

        if(!is.matrix(tyw)) {
            tyw <- matrix(tyw)
            tww <- array(tww,dim = c(ncol(W),ncol(W),n))
        }

        nc <- ncol(tww[,,1])

        coef.mat <- matrix(cv.maxPenalty,ncol(W),n)
        epsilon <- 1.0/n
        ridge <- double(n)
        ridge.lc <- double(n)
        doridge <- !logical(n)

        ridger <- function(i) {
            doridge[i] <<- FALSE
            ridge.lc[i] <- ridge[i]*tyw[1,i][1]/NZD(tww[,,i][1,1])
            tryCatch(chol2inv(chol(tww[,,i]+diag(rep(ridge[i],nc))))%*%tyw[,i],
                     error = function(e){
                         ridge[i] <<- ridge[i]+epsilon
                         doridge[i] <<- TRUE
                         return(rep(cv.maxPenalty,nc))
                     })
        }

        while(any(doridge)){
            iloo <- (1:n)[doridge]
            coef.mat[,iloo] <- sapply(iloo, ridger)
        }

        m <- sapply(1:n, function(i) {
            W[i,, drop = FALSE] %*% coef.mat[,i]
        })

    }

    ## Check for finite m/objective function, penalize if not

    if(all(is.finite(m))) {
        return(mean(check.function(y-m,tau=tau)))
    } else {
        return(cv.maxPenalty)
    }

}

## Data-driven choice of delta and bandwidths for x. The function
## returns a list with elements delta (scalar) and bws (vector if
## multivariate predictors are used).

npqw.optim <- function(x,
                       y,
                       tau=0.5,
                       nmulti=1,
                       sd.pilot,
                       ckertype="gaussian",
                       regtype="lc") {

    if(missing(x)) stop("must provide x")
    if(missing(y)) stop("must provide y")
    if(missing(sd.pilot)) stop("must provide sd.pilot")
    if(tau < 0.0001 || tau > 0.9999) stop("tau must lie in [0.0001,0.9999]")

    x <- as.data.frame(x)

    n <- length(y)
    x.col.numeric <- sapply(1:NCOL(x),function(i){is.numeric(x[,i])})

    fv <- cv.maxPenalty <- sqrt(.Machine$double.xmax)

    for(iMulti in 1:nmulti) {

        optim.return <- list()
        optim.return$value <- cv.maxPenalty

        while(optim.return$value==cv.maxPenalty) {

            ## We need to feed the initial values for optimization to
            ## optim along with constraints. As we may have
            ## mixed-data, we check which variables are numeric and
            ## use initial values guided by rules-of thumb. When
            ## nmulti > 1 we add noise to these starting values.

            init <- numeric(1+NCOL(x))
            lower <- init
            upper <- init

            if(iMulti==1) {
                init[1] <- tau
                for(i in 1:ncol(x)) {
                    if(x.col.numeric[i]) init[i+1] <- sd(x[,i])*n^{-.2}
                    if(!x.col.numeric[i]) init[i+1] <- 0.25
                }
            } else {
                init[1] <- runif(1)
                for(i in 1:ncol(x)) {
                    if(x.col.numeric[i]) init[i+1] <- runif(1,0.1,2.5)*sd(x[,i])*n^{-.2}
                    if(!x.col.numeric[i]) init[i+1] <- runif(1,0.1,0.5)
                }
            }

            lower[1] <- .Machine$double.eps
            upper[1] <- 1-.Machine$double.eps
            for(i in 1:ncol(x)) {
                if(x.col.numeric[i]) {
                    lower[i+1] <- 0.001*sd(x[,i])*n^{-.2} 
                    lower[i+1] <- 0.1*sd(x[,i])*n^{-.2} ## Changed to require smoother fits HLB
                    upper[i+1] <- .Machine$double.xmax
                }
                if(!x.col.numeric[i]) {
                    lower[i+1] <- 0
                    upper[i+1] <- 1
                }
            }

            optim.return <- optim(init,
                                  npqw.estimator.loo,
                                  x=x,
                                  y=y,
                                  tau=tau,
                                  ckertype=ckertype,
                                  regtype=regtype,
                                  sd.pilot=sd.pilot,
                                  lower=lower,
                                  upper=upper,
                                  method="L-BFGS-B",
                                  control = list(maxit = 10000))

        }

        print(paste("tau = ",tau,
                    ": restart ",iMulti,
                    " of ",nmulti,
                    ": fval = ",formatC(optim.return$value,format="f",digits=4),sep=""))

        if(optim.return$value < fv) {
            fv <- optim.return$value
            optim.out <- optim.return
        }

    }

    return(list(delta=optim.out$par[1],bws=optim.out$par[-1]))

}

