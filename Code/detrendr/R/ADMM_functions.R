
#' Proximal mapping of quantile trend filtering objective
#'
#' \code{prox_f} computes the proximal mapping of the quantile trend filtering objective
#'
#' @param x0 input initial value of c(theta, eta)
#' @param y observed data
#' @param tau quantile level
#' @param lambda penalty parameter for differences
#' @param step_size step parameter of prox
prox_f <- function(x0, step_size, lambda, y, D, M, tau){
  n <- length(y)
  x <- x0
  x[1:n] <- prox_f1(x0[1:n], y, tau, step_size)
  x[(n+1):length(x)] <- prox_f2(x0[(n+1):length(x)], lambda, step_size)
  return(x)
}

#' Proximal mapping of indicator eta = D*theta
#'
#' \code{prox_g} Computes the proximal mapping of the constraint (projection)
#'
#' @param x0 input initial value of c(theta, eta)
#' @param D differencing matrix
#' @param M cholesky factorization of (crossprod(D) + I)
prox_g <- function(x0, lambda, y, D, M, tau){
  n <- ncol(D)
  proj_x <- project_V_R(x0[1:n], x0[(n+1):length(x0)], D, M)
  return(c(proj_x$theta, proj_x$eta))
}

#' Proximal ADMM algorithm
#'
#' \code{ADMM} Runs the proximal admm algorithm for minimizing
#' f1(x) + f2(z) subject to z = Dx.
#'
#' @param x0 initial value of primal variable (c(x0, z0))
#' @param y0 initial value of dual variable
#' @param prox_f function for computing the proximal mapping of the objective
#' @param prox_g function for computing the proximal mapping of the constraint
#' @param step_size step parameter for prox_f
#' @param max_iter maximum number of iterations
#' @param eps_abs absolute epsilon for stopping criteria
#' @param eps_rel relative epsilon for stopping criteria
ADMM <- function(x0, y0, prox_f, prox_g, step_size,
                 max_iter, eps_abs, eps_rel,
                 ...){


  n <- length(y)
  dual_norm <- double(max_iter)
  primal_norm <- double(max_iter)
  eps_pri <- double(max_iter)
  eps_dual <- double(max_iter)
  z0 <- prox_g(x0 + y0, ...)
  x_k <- prox_f(z0-y0, step_size, ...)
  z_k <- prox_g(x_k + y0, ...)
  y_k <- (y0 + x_k - z_k)

  sqrt_n <- sqrt(n)

  for (iter in 1:max_iter){
    x_k <- prox_f(z_k-y_k, step_size, ...)
    dual_resid <- -z_k
    z_k <- prox_g(x_k + y_k, ...)
    dual_resid <- -(dual_resid + z_k)/step_size
    y_k <- (y_k + x_k - z_k)
    primal_resid <- x_k - z_k
    dual_norm[iter] <- norm(dual_resid, type = "2")
    primal_norm[iter] <- norm(primal_resid, type = "2")

    eps_pri[iter] <- sqrt_n*eps_abs +
      eps_rel*max(norm(x_k, type = "2"),
                  norm(z_k, type = "2"))
    eps_dual[iter] <- sqrt_n*eps_abs + eps_rel*norm(y_k/step_size, type = "2")

    if(primal_norm[iter] < eps_pri[iter] &
       dual_norm[iter] < eps_dual[iter]) {
      break
    }

  }
  return(list(x = x_k, y = y_k,
              dual_norm = dual_norm,
              primal_norm = primal_norm,
              eps_pri = eps_pri, eps_dual=eps_dual))
}



#' Linearized ADMM algorithm
#'
#' \code{ADMM} Runs the Linearized ADMM algorithm for minimizing
#' f(x) + g(Ax)
#'
#' @param x0 initial value of primal variable
#' @param y0 initial value of dual variable
#' @param prox_f function for computing the proximal mapping of f
#' @param prox_g function for computing the proximal mapping of g
#' @param step_f step parameter for prox_f with 0 < step_f < step_g/||A||^2
#' @param step_g step parameter for prox_g
#' @param max_iter maximum number of iterations
LADMM <- function(x0, y0, A, prox_f, prox_g, step_f, step_g,
                 max_iter, ...){
  z0 <- A%*%x0

  x_k <- prox_f(x0 - (step_f/step_g) * crossprod(A, A%*%x0 - z0 + y0),
                step_f, ...)
  z_k <- prox_g(A%*%x_k + y0, step_g, ...)
  y_k <- (y0 + A%*%x_k - z_k)

  for (iter in 1:max_iter){
    x_k <- prox_f(x_k - (step_f/step_g) * crossprod(A, A%*%x_k - z_k + y_k),
                  step_f, ...)
    z_k <- prox_g(A%*%x_k + y_k, step_g, ...)
    y_k <- (y_k + A%*%x_k - z_k)
  }

  return(list(x = x_k, y = y_k))
}


