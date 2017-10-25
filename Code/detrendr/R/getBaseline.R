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
  k <- 4
  n <- length(y)
  theta <- y
  D <- get_Dk(n, k)
  eta <- matrix(D%*%theta)
  M <- Diagonal(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)
  lambda <- lambda0*n^(k-1)/factorial(k-1)
  theta <- warmStart(y, k, lambda0, tau, 5)
  multi_step <- spingarn_multi_step(theta, eta, y, D, cholM,
                                    lambda, tau, (max(y)-min(y))/5, maxiter, k)
  plot(log(multi_step[[3]]), type="l")
  theta <- multi_step[[1]]
  theta_last <- prox_f1(theta, y, tau)
  theta_last[which(is.na(y0))] <- NA
  plot(y, type="l")
  lines(theta_last, col="red")
  return(theta_last)
}
