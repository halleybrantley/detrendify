#' @title 
#' Get baseline using trend filtering
#' @description 
#' Get smooth baseline of data. Prints data and baseline estimate every 5000
#' iterarions. Data may contain missing values. Missing values are replaced using 
#' next observation carry backward. 
#' @param y vector of numbers to fit baseline
#' @param lambda regularization parameter, larger values produce 
#' smoother baselines
#' @param tau quantile level
#' @param maxiter number of iterations to use in algorithm
#' @examples 
#' 
#' @export

getBaseline <- function(y, lambda0 = 1e-10, maxiter = 20000, tau=0.05){
  y0 <- y
  if(is.na(y[length(y)])) {y[length(y)] <- y[max(which(!is.na(y)))]}
  y <- zoo::na.locf(y, fromLast = TRUE)
  k <- 3
  n <- length(y)
  D <- get_Dk(n, k)
  theta <- warmStart(y, k, lambda0, tau, 5)
  eta <- matrix(D%*%theta)
  M <- Diagonal(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)

  # theta <- y
  # eta <- matrix(D%*%theta)
  
  lambda <- lambda0*n^2/(1000^2)
  step <- .2
  multi_step <- spingarn_multi_step(theta, eta, y, D, cholM,
                                     lambda, tau, step, 10000, k)
  

  theta <- multi_step[[1]]
  eta <- multi_step[[2]]
  theta_last <- prox_f1(theta, y, tau, step)
  eta_last <- prox_f2(eta, lambda, step)
  
  theta_last[which(is.na(y0))] <- NA
  plot(y, type="l")
  lines(theta_last, col="red")
  
  return(theta_last)
}
