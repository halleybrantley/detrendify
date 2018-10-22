
# FASTA functions

# Differentiable Objective
f <- function(x, y, tau, D, nu, lambda){
  n <- length(y)
  theta <- x[1:n]
  e <- y - theta
  rho <-  (log(cosh(e*25*pi))/2)/(25*pi) + (tau - .5)*e
  eta <- x[(n+1):length(x)]
  penalty <- sum((eta-D%*%theta)^2)/(2*nu)
  return(sum(rho) + penalty)
}

# Non-differentiable Objective
g <- function(x, y, tau, D, nu, lambda){
  n <- length(y)
  eta <- x[(n+1):length(x)]
  return(lambda*sum(abs(eta)))
}

# Gradient of differentiable part
gradf <- function(x, y, tau, D, nu, lambda){
  n <- length(y)
  theta <- x[1:n]
  eta <- x[(n+1):length(x)]
  gf_eta <- (eta - D%*%theta)/nu
  gf_theta <- -(tanh((y-theta)*25*pi)/2 - 1/2 + tau) -
    crossprod(D, gf_eta)
  return(rbind(gf_theta, gf_eta))
}

# Prox operator of non-differentiable part
proxg <- function(x, step_size, y, tau, D, nu, lambda){
  n <- ncol(D)
  prox <- x
  eta <- x[(n+1):length(x)]
  eta_p <- double(length(eta))
  eta_p[eta > lambda*step_size] <- eta[eta > lambda*step_size] - lambda*step_size
  eta_p[eta < -lambda*step_size] <- eta[eta < -lambda*step_size] + lambda*step_size
  prox[(n+1):length(x)] <- eta_p
  return(prox)
}

