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

getBaseline <- function(y, lambda0 = 1e-10, tau=0.05){
  y0 <- y
  if(is.na(y[length(y)])) {y[length(y)] <- y[max(which(!is.na(y)))]}
  y <- zoo::na.locf(y, fromLast = TRUE)
  k <- 3
  n <- length(y)
  D <- get_Dk(n, k)
  lambda <- lambda0
  theta <- warmStart(y, k, lambda, 5/n, tau, 5, n*600/5)
  
  eta <- matrix(D%*%theta)
  M <- Diagonal(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)

  multi_step <- spingarn_multi_step(theta, eta, y, D, 
                                    cholM, lambda, tau, 1/n, n*50, k)
  
  # plot(diff(multi_step[[3]][20000:maxiter]), type="l")
  # plot(multi_step[[4]][100:maxiter], type="l")
  
  theta_last <- prox_f1( multi_step[[1]], y, tau, step)
  eta_last <- prox_f2(multi_step[[2]], lambda, step)
  
  theta_last[which(is.na(y0))] <- NA
  plot(y, type="l")
  lines(theta_last, col="blue")
  lines(theta, col="red")
  return(theta_last)
}
