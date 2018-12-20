#' Get quantile trend using Spingarn's Method
#'
#' \code{get_trend_spingarn} Returns the quantile trend vector
#'
#' @param y observed data, should be equally spaced, may contain NA
#' @param tau quantile level at which to evaluate trend
#' @param lambda penalty paramter controlling smoothness
#' @param k order of differencing
#' @export
get_trend_spingarn <- function(y, tau, lambda, k, numIter = 10000){
  if (length(tau)!=1){
    stop("Can only evaluate single quantile")
  }
  
  theta0 <- y
  D <- get_Dk(length(y), k)
  eta0 <- as.numeric(D%*%theta0)
  M <- Matrix::chol(Diagonal(n) + Matrix::crossprod(D))
  spign_const <- spingarn_multi_step(theta0, eta0, y, D, M, lambda, 
                                     tau, 1, numberIter=numIter, k)
  return(spign_const[["theta"]])
}
  