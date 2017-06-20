#' Quantile Proximal mapping
#'
#' \code{prox_quantile_R} computes the proximal mapping of the check function.
#'
#' @param w input
#' @param tau quantile parameter
#' @param alpha scale parameter
#' @examples
#' set.seed(12345)
#' n <- 1e3
#' w <- seq(-3, 3, length.out=n)
#' tau <- 0.5
#' alpha <- 2
#' prox_out <- prox_quantile(w, tau, alpha)
#' plot(w, prox_out, type='l', main=expression(paste(tau," = ")))
#'
#' tau <- 0.05
#' alpha <- 2
#' prox_out <- prox_quantile(w, tau, alpha)
#' plot(w, prox_out, type='l', main=expression(paste(tau," = ")))
prox_quantile_R <- function(w, tau, alpha) {
  n <- length(w)
  prox_out <- double(n)
  threshold1 <- tau*alpha
  threshold2 <- -(1 - tau)*alpha
  ix1 <- which(w > threshold1)
  ix2 <- which(w < threshold2)
  prox_out[ix1] <- w[ix1] - threshold1
  prox_out[ix2] <- w[ix2] - threshold2
  return(matrix(prox_out, ncol=1))
}

#' Proximal mapping of f_1
#'
#' \code{prox_f1} computes the proximal mapping of the average quantile loss
#'
#' @param theta input
#' @param y response
#' @param tau quantile parameter
#' @param t step-size
#' set.seed(12345)
prox_f1_R <- function(theta, y, tau = 0.05, t = 1) {
  n <- length(theta)
  w <- y - theta
  alpha <- t/n
  return(y - prox_quantile(w, tau, alpha))
}

#' Proximal mapping of f_2
#'
#' \code{prox_f2} computes the proximal mapping of the average quantile loss
#'
#' @param eta input
#' @param lambda regularization parameter
#' @param t step-size
#' @examples
#' set.seed(12345)
#' n <- 1e3
#' eta <- seq(-3, 3, length.out=n)
#' lambda <- 1
#' prox_out <- prox_f2(eta, lambda)
#' plot(eta, prox_out, type = 'l')
#' abline(0,1)
prox_f2_R <- function(eta, lambda, t = 1) {
  return(prox_quantile(eta, tau = 0.5, alpha = 2*t*lambda))
}

#' Proximal mapping
#'
#' \code{prox} computes the block separable proximal mapping.
#'
#' @param theta input
#' @param eta input
#' @param y response
#' @param lambda regularization parameter
#' @param tau quantile parameter
#' @param t step-size
#' @examples
#' set.seed(12345)
prox_R <- function(theta, eta, y, lambda, tau = 0.05, t=1) {
  return(list(theta=prox_f1(theta, y, tau, t), 
              eta = prox_f2(eta, lambda, t)))
}

#' First order difference matrix
#'
#' \code{get_D1_R} computes the discrete derivative matrix.
#'
#' @param n length of input
#' @examples
#' n <- 5
#' D1 <- get_D1(n)
get_D1_R <- function(n) {
  D1 <- Matrix(data=0, nrow=n-1, ncol=n, sparse = TRUE)
  for (i in 1:(n-1)) {
    D1[i, i] <- 1
    D1[i, i+1] <- -1
  }
  return(D1)
}

#' kth order difference matrix
#'
#' \code{get_Dk_R} computes the discrete kth derivative matrix.
#'
#' @param n length of input
#' @param k order of the derivative
#' @examples
#' n <- 10
#' k <- 3
#' D <- get_Dk(n, k)
get_Dk_R <- function(n, k) {
  D <- get_D1_R(n)
  for (i in 2:k) {
    D <- get_D1_R(n - i + 1) %*% D
  }
  return(D)
}

#' Project onto subspace
#'
#' \code{project_V} projects (theta, eta) onto the subspace eta = D%*%theta
#'
#' @param theta first input
#' @param eta second input
#' @param D differencing matrix
#' @examples
#' set.seed(12345)
#' k <- 3
#' n <- 1e2
#' D <- get_Dk(n, k)
#' theta <- rnorm(n)
#' eta <- D %*% theta
#' proj <- project_V(theta, eta, D)
#'
#' m <- 1e3
#' angles <- double(m)
#' theta2 <- matrix(rnorm(n*m), n, m)
#' eta2 <- D%*%theta2
#' angles <- matrix(t(theta - proj$theta)%*%(theta2 - matrix(proj$theta)%*%matrix(1,1,m)))
#' angles <- angles + matrix(t(eta - proj$eta)%*%(eta2 - matrix(proj$eta)%*%matrix(1,1,m)))
#' summary(angles)
project_V_R <- function(theta, eta, D, M) {
  theta <- Matrix::solve(M, theta + Matrix::t(D)%*%eta, sparse=TRUE)
  eta <- D%*%theta
  eta <- matrix(eta, ncol=1)
  theta <- matrix(theta, ncol=1)
  return(list(theta=theta, eta=eta))
}

#' One step of Spingarn's algorithm
#'
#' \code{spingarn_one_step_R} computes a single round of Spingarn updates
#' using R code.
#'
#' @param theta input 1
#' @param eta input 2
#' @param y response
#' @param D differencing matrix
#' @param lambda regularization parameter
#' @param tau quantile parameter
#' @param t step-size
#' @export
#' @examples
#' set.seed(12345)
#' n <- 1e2
#' x <- seq(1/n, 1, length.out=n)
#' f <- 2*(x + 2)^2 + 3*cos(3*pi*x)
#' tau <- 1e4
#' g <-100*exp(-tau*(x-0.5)^2)
#' y <- f + g + rnorm(n)
#' k <- 3
#' D <- get_Dk_R(n, k)
#' M <- diag(n) + crossprod(D)
#' lambda <- 1
#' tau <- 0.01
#' t <- 1
#' theta <- y
#' eta <- matrix(D %*% theta)
#'
#' max_iter <- 1e2
#' THx <- matrix(NA, n, max_iter)
#' for (iter in 1:max_iter) {
#'    one_step <- spingarn_one_step_R(theta, eta, y, D, M, lambda, tau=tau, t=t)
#'    theta <- one_step$theta
#'    eta <- one_step$eta
#'    THx[,iter] <- prox_f1(theta, y, tau)
#' }
#' theta_last <- prox_f1(theta, y, tau, t=t)
#' rerr <- double(max_iter-1)
#' for (i in 1:(max_iter-1)) {
#'   rerr[i] <- norm(as.matrix(THx[,i+1] - THx[,i]),'f')/(1 + norm(THx[,i,drop=FALSE], 'f'))
#' }
#' plot(x,f,type='l',col='blue', ylim=c(min(y),max(y)), lwd=3)
#' points(x,y,pch=16)
#' lines(x,theta_last,col='red', lwd=3)
spingarn_one_step_R <- function(theta, eta, y, D, M, lambda, tau=0.05, t=1) {
  theta_old <- theta
  eta_old <- eta
  prox_sol <- prox_R(theta, eta, y, lambda, tau, t)
  theta <- prox_sol$theta
  eta <- prox_sol$eta
  proj_sol <- project_V_R(2*theta - theta_old, 2*eta - eta_old, D, M)
  theta <- theta_old + 1.9*(proj_sol$theta - theta)
  eta <- eta_old + 1.9*(proj_sol$eta - eta)
  return(list(theta=theta, eta=matrix(eta)))
}

#' Multiple step of Spingarn's algorithm
#'
#' \code{spingarn_one_step} computes a single round of Spingarn updates.
#'
#' @param theta input 1
#' @param eta input 2
#' @param y response
#' @param D differencing matrix
#' @param lambda regularization parameter
#' @param tau quantile parameter
#' @param t step-size
#' @export
spingarn_multi_step_R <- function(theta, eta, y, D, M, lambda, 
                                  tau=0.05, t=1, numberIter=1){
  for (iter in 1:numberIter) {
    one_step <- spingarn_one_step_R(theta, eta, y, D, M, lambda, tau, step)
    theta <- one_step[[1]]
    eta <- one_step[[2]]
  }
  return(list(theta=theta, eta=matrix(eta)))
}
