generate_peaks_design <- function(n){
  x <- seq(0.5, n, 1)/n
  df <- sample((n/50):(n/25), 1)
  splineBasis <- ns(x, df=df)
  theta <- rexp(ncol(splineBasis), 1)
  baseline <- splineBasis%*%theta
  # plot(baseline[1:300], type="l")
  
  numberOfPeaks <- sample((n/100):(n/50), 1)
  peakCenters <- runif(numberOfPeaks)*.8+.1
  peakArea <- runif(numberOfPeaks)*20/n
  peakWidths <- runif(numberOfPeaks)*4/n + 1/n
  
  peaks <- rep(0, length(baseline))
  if (numberOfPeaks > 0){
    for (i in 1:numberOfPeaks){
      peaks <- peaks + 
        peakArea[i]*dnorm(x, mean=peakCenters[i], sd = peakWidths[i])
    }
  }
  # plot(peaks[1:300], type="l")
  
  noise <- 0.25*rnorm(n)
  y <- peaks + baseline + noise
  # plot(y~x, type="l", col="grey")
  # lines((peaks+baseline)~x)
  # lines(baseline~x, col="red")
  
  df <- data.frame(y=y, x=x, baseline=baseline, peaks = peaks)
  return(df)
}

generate_gaus <- function(n){
  x <- seq(0.5, n, 1)/n
  f <- sin(2*pi*x) 
  y <- f + ((1+x^2)/4)*rnorm(n)
  df <- data.frame(y=y, x=x, f=f)
  return(df)
}

generate_t <- function(n){
  x <- seq(0.5, n, 1)/n
  f <- sin(2*pi*x) 
  y <- f + rt(n, df=5)
  df <- data.frame(y=y, x=x, f=f)
  return(df)
}

generate_exp <- function(n){
  x <- seq(0.5, n, 1)/n
  f <- sin(2*pi*x) 
  y <- f + rexp(n, 1)
  df <- data.frame(y=y, x=x, f=f)
  return(df)
}


generate_shapebeta <- function(n){
  x <- seq(0.5, n, 1)/n
  f <- sin(2*pi*x) 
  y <- f + rbeta(n, 1, 11-10*x)
  df <- data.frame(y=y, x=x, f=f)
  return(df)
}

generate_mixednorm <- function(n){
  x <- seq(0.5, n, 1)/n
  f <- sin(2*pi*x) 
  e1 <- rnorm(n, 1, sqrt(.25))
  e2 <- rnorm(n, -1, sqrt(.25))
  draw <- rbinom(n, 1, x)
  y <- f + draw*e1 + (1-draw)*e2
  df <- data.frame(y=y, x=x, f=f)
  return(df)
}

trueQuantile <- function(simDesign, x, tau){
  Q <- matrix(NA, length(x), length(tau))
  for (i in 1:length(tau)){
    if (simDesign == "gaus"){
      Q[,i] <- sin(2*pi*x) + ((1+x^2)/4)*qnorm(tau[i])
    } else if (simDesign == "shapebeta") {
      Q[,i] <- sin(2*pi*x) + qbeta(tau[i], 1, 11-10*x)
    } else if (simDesign == "mixednorm") {
      cdf <- function(e, x, tau){
        x*pnorm(e, 1, sqrt(.25)) + (1-x)*pnorm(e, -1, sqrt(.25)) - tau
      }
      q <- sapply(x, function(x) uniroot(cdf, c(-5, 5), x=x, tau = tau[i])$root)
      Q[,i] <- sin(2*pi*x) + q
    }
  }
  return(Q)
}
