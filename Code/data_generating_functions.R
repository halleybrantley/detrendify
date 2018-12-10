generate_peaks_design <- function(n){
  x <- seq(0.5, n, 1)/n
  df <- sample(5:20, 1)
  splineBasis <- ns(x, df=df)
  theta <- rexp(ncol(splineBasis), .5)
  baseline <- splineBasis%*%theta
  
  plot(splineBasis[,1], type="l")
  for(i in 1:ncol(splineBasis)){
    lines(splineBasis[,i])
  }
  plot(baseline, type="l")
  
  numberOfPeaks <- sample(5:15, 1)
  peakCenters <- runif(numberOfPeaks)*.8+.1
  peakArea <- runif(numberOfPeaks)*.05
  peakWidths <- runif(numberOfPeaks)*.01 + 2/n
  
  peaks <- rep(0, length(baseline))
  if (numberOfPeaks > 0){
    for (i in 1:numberOfPeaks){
      peaks <- peaks + 
        peakArea[i]*dnorm(x, mean=peakCenters[i], sd = peakWidths[i])
    }
  }
  noise <- 0.25*rnorm(n)
  y <- peaks + baseline + noise
  plot(y~x, type="l", col="grey")
  lines((peaks+baseline)~x)
  lines(baseline~x, col="red")
  
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