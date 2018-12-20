

trend_spingarn <- function(y, tau, lambda, k){
  n <- length(y)
  D <- get_Dk(n, k)
  theta <- rep(quantile(y, tau), n)
  eta <- D%*%theta
  M <- Matrix::Diagonal(n) + Matrix::crossprod(D)
  cholM <- Matrix::chol(M)

  # result <- spingarn_multi_step(theta, eta, y,  D, cholM,
  #                     lambda, tau, step = 1,
  #                     numberIter=1, k, rel_tol = 0.0001)
}