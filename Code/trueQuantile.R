trueQuantile <- function(simDesign, x, tau){
  Q <- matrix(NA, length(x), length(tau))
  for (i in 1:length(tau)){
    if (simDesign == "gaus"){
      Q[,i] <- sin(2*pi*x) + ((1+x^2)/4)*qnorm(tau[i])
    } else if (simDesign == "shapebeta") {
      Q[,i] <- sin(2*pi*x) + qbeta(tau[i], 1, 11-10*x)
    } else if (simDesign == "mixednorm") {
      cdf <- function(e, x, tau){
        x*pnorm(e, -1, sqrt(.25)) + (1-x)*pnorm(e, 1, sqrt(.25)) - tau
      }
      q <- sapply(x, function(x) uniroot(cdf, c(-5, 5), x=x, tau = tau[i])$root)
      Q[,i] <- sin(2*pi*x) + q
    }
  }
  return(Q)
}