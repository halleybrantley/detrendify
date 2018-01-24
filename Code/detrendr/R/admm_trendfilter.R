

censor <- function(y){
  z <- y
  z[z<0] <- 0
  return(z)
}

admm_trendfilter <- function(y, tau, lambda, k, rho, maxiter){
  n <- length(y)
  D <- get_Dk(length(y), k)
  m <- nrow(D)
  A <- cbind(D, -D, -diag(m), diag(m))
  b <- as.numeric(D%*%y)
  
  theta <- rep(mean(y), n)   
  eta <- rep(0, m)
  r <- y-theta
  r_plus <- censor(r)
  r_minus <- censor(-r)
  eta_plus <- censor(eta)
  eta_minus <- censor(-eta)
  
  q <- c(rep(tau, n), rep(1-tau, n), rep(lambda, 2*m))
  M <- rbind(cbind(rho*diag(ncol(A)), Matrix::t(A)), cbind(A, 0*diag(nrow(A)))) 

  x <- c(r_plus, r_minus, eta_plus, eta_minus)  
  z <- x
  u <- rep(0, length(x))
  
  for (i in 1:maxiter){
    rhs <- -c(q-rho*(z-u), -b)
    sol <- Matrix::solve(M, rhs)
    x <- sol[1:(2*n+2*m)]
    z <- censor(x+u)
    u <- u + x - z
  }
  
  theta <- y - x[1:n] + x[(n+1):(2*n)]
  return(theta)
}