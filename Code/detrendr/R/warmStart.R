#' @title 
#' Warm start for Spingarn algorithm 
#' @description 
#' Get approximation of initial values by aggregating data and fitting on 
#' smaller dataset
#' @param y original data to fit baseline
#' @param k order of differencing matrix
#' @param lambda regularization parameter
#' @param tau quantile level
#' @param reduction number of observations to aggregate
#' @export
warmStart <- function(y, k, lambda, step, tau, reduction, maxiter){
  require(dplyr)
  x <- seq(1, length(y), 1)
  df <- data.frame(y=y, x=x)
  df$xBreaks <- cut(df$x, seq(0, length(x)+reduction, reduction))
  dfAgg <- ddply(df, .(xBreaks), summarize, 
                 yMean = quantile(y, probs=get("tau")))
  
  y2 <- dfAgg$yMean
  theta <- y2
  n <- length(y2)
  D <- get_Dk(n, k)
  M <- Diagonal(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)
  eta <- matrix(D%*%theta)
  
  multi_step <- spingarn_multi_step(theta, eta, y2, D, cholM,
                                    lambda, tau, step, 
                                    maxiter, k)
  #print(paste("converged in", multi_step[[6]], "iterations."))
  #plot(multi_step[[3]], type="l")
  
  dfAgg$theta <- prox_f1(multi_step[[1]], y2, tau)
  plot(y2, type="l")
  lines(dfAgg$theta, col="red")

  dfAll <- merge(df, dfAgg, by = "xBreaks")
  dfAll <- dfAll[order(dfAll$x), ]
  return(dfAll$theta)
}
