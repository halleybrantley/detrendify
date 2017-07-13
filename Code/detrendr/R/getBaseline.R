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
getBaseline <- function(y, lambda0 = 1e-7, maxiter = 20000){
  y0 <- y
  if(is.na(y[length(y)])) {y[length(y)] <- y[max(which(!is.na(y)))]}
  y <- zoo::na.locf(y, fromLast = TRUE)
  tau <- .02
  k <- 3
  n <- length(y)
  theta <- y
  D <- get_Dk(n, k)
  eta <- matrix(D%*%theta)
  M <- Diagonal(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)
  lambda <- lambda0*n^(k-1)/factorial(k-1)
  theta <- warmStart(y, k, lambda0, tau, 10)
  
  pb <- txtProgressBar(min = 1, max = (maxiter/5000))
  for (i in 1:(maxiter/5000)){
    multi_step <- spingarn_multi_step(theta, eta, y, D, cholM,
                                    lambda, tau, 1, 5000, k)
    setTxtProgressBar(pb, i)
    theta <- multi_step[[1]]
    eta <- multi_step[[2]]
    theta_last <- prox_f1(theta, y, tau)
    plot(y, type="l")
    lines(theta_last, col="red")
  }
  close(pb)
  theta <- multi_step[[1]]
  theta_last <- prox_f1(theta, y, tau)
  theta_last[which(is.na(y0))] <- NA
  return(theta_last)
}
