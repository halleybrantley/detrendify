library(splines)
library(fields)
n <- 21000
tau <- c(0.01, .05, .1)
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

numberOfPeaks <- sample(10:25, 1)
peakCenters <- runif(numberOfPeaks)*.8+.1
peakArea <- runif(numberOfPeaks)*.05
peakWidths <- runif(numberOfPeaks)*.01 + 100/n

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
lines((baseline+qnorm(.1, sd = .25))~x, col="red")

df <- data.frame(y=y, x=x, baseline=baseline, peaks = peaks)


trend <- matrix(NA, n, length(tau))
for (j in 1:length(tau)){
  fit_qsreg <- qsreg(df$x, df$y, maxit.cv = 50, alpha=tau[j], 
                     hmin = -30)
  trend[,j] <- predict(fit_qsreg)    
}

plot(y~x, type="l", col="grey")
lines((baseline+qnorm(.1, sd = .25))~x, col="red")
for (i in 1:length(tau)){
  lines(trend[,i]~df$x, col="blue")
}


lam <- lambda_eBIC(df$y, tau, 3,  gamma = 1,
                   lambdaSeq = n^seq(1, 1.4, length.out=15),
                   plot_lambda = FALSE, single_lambda = FALSE)
trend2 <- gurobi_trend(df$y, tau, lam$lambda, k=3)

for (i in 1:length(tau)){
  lines(trend2[,i]~df$x, col="darkgreen")
}

 mean((trend[,1] - (df$baseline+qnorm(tau[1], sd = 0.25)))^2)
 mean((trend[,2] - (df$baseline+qnorm(tau[2], sd = 0.25)))^2)
 mean((trend[,3] - (df$baseline+qnorm(tau[3], sd = 0.25)))^2)
 
 mean((trend2[,1] - (df$baseline+qnorm(tau[1], sd = 0.25)))^2)
 mean((trend2[,2] - (df$baseline+qnorm(tau[2], sd = 0.25)))^2)
 mean((trend2[,3] - (df$baseline+qnorm(tau[3], sd = 0.25)))^2)
 