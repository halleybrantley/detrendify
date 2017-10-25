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

getBaseline <- function(y, lambda0 = 0.1, maxiter = 20000, tau=0.05){
  y0 <- y
  if(is.na(y[length(y)])) {y[length(y)] <- y[max(which(!is.na(y)))]}
  y <- zoo::na.locf(y, fromLast = TRUE)
  k <- 4
  n <- length(y)
  theta <- y
  D <- get_Dk(n, k)
  eta <- matrix(D%*%theta)
  M <- Diagonal(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)
  lambda <- lambda0*n^(k-1)/factorial(k-1)
  multi_step <- spingarn_multi_step(theta, eta, y, D, cholM,
                                  lambda, tau, 1, 50000, k)

  theta <- multi_step[[1]]
  theta_last <- prox_f1(theta, y, tau)
  plot(y, type="l")
  lines(theta_last, col="red")

  return(theta_last)
}
