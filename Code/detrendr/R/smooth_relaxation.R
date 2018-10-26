
approx_checkLoss <- function(x, tau){
  x <- as.numeric(x)
  f <- (log(cosh(x*1e4*pi))/(2*1e4*pi) + (tau - .5)*x)
  f[is.infinite(f) & x > 0] <- tau*x[is.infinite(f) & x > 0]
  f[is.infinite(f) & x < 0] <- (tau-1)*x[is.infinite(f) & x < 0]
  return(f)
}

hess_checkLoss <- function(x){
  x <- as.numeric(x)
  return(1e4*pi*cosh(1e4*pi*x)^(-2)/2)
}

grad_checkLoss <- function(x, tau){
  x <- as.numeric(x)
  return(tanh(x*1e4*pi)/2 - 1/2 + tau)
}

obj <- function(phi, y, tau, D, lambda){
  return(
    sum(approx_checkLoss(phi, tau),
    lambda*approx_checkLoss(D%*%(y-phi), 0.5))
    )
}

grad_obj <- function(phi, y, tau, D, lambda){
  as.numeric(grad_checkLoss(phi, tau) -
    lambda*Matrix::crossprod(D,grad_checkLoss(D%*%(y-phi), 0.5)))
}

hess_obj <- function(phi, y, tau, D, lambda){
  diag(hess_checkLoss(phi)) +
    lambda*Matrix::crossprod(D, diag(hess_checkLoss(D%*%(y-phi))))%*%D
}

