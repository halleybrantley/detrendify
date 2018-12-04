#' Consensus ADMM for overlapping windows
#'
#' \code{consensus ADMM} Uses ADMM to smooth overlapping windows
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param lambda smoothing penalty parameter
#' @param k order of differencing
#' @param rho step size for ADMM
#' @param window_size size of windows to use
#' @param overlap integer length of overlap between windows
#' @param max_iter Maximum number of iterations
#' @examples
#' n <- 100
#' x <- seq(1, n, 1)
#' y <- sin(x*2*pi/n) + rnorm(n, 0, .4)
#' lambda <- 10
#' k <- 3
#' tau <- c(0.05)
#' overlap <- 20
#' rho <- 1
#' result <- consensus_ADMM(y, tau, lambda, k, rho, overlap, 300)
#' y_n <- length(y)
#' window_size <- floor((y_n+overlap)/2)
#' plot(result$phiBar~x, type="l")
#' points(result$phi1~x[1:window_size])
#' points(result$phi2~x[(y_n-window_size+1):y_n], col="blue")
#' plot(y~x)
#' lines(result$theta)
#' plot(result$primal_norm)
#' plot(result$dual_norm)
#' @export
consensus_ADMM <- function(y, tau, lambda, k, rho, window_size,
                           overlap, max_iter, eps = 0.01, update=10){
  y_n <- length(y)
  tau <- sort(tau)
  nT <- length(tau)
  D <- get_Dk(window_size, k)
  n_windows <- ceiling((y_n-overlap)/(window_size-overlap))
  windows <- matrix(FALSE, y_n, n_windows)
  y_list <- list()
  w_list <- list()
  phi_list <- list()
  
  # Initial values

  for (i in 1:n_windows){
    start_I <- 1+(window_size-overlap)*(i-1)
    end_I <- min((window_size + (window_size-overlap)*(i-1)), y_n)
    windows[start_I:end_I,i] <- TRUE
    y_list[[i]] <- y[start_I:end_I]
    w_list[[i]] <- numeric(length(y_list[[i]])*length(tau))
    phi_list[[i]] <- matrix(0,length(y_list[[i]]),length(tau))
  }

  model_list <- lapply(y_list, create_model, tau, lambda, D, rho)
  
  overlapInd <- rowSums(windows) > 1
  phiBar_list <- update_consensus(phi_list, windows, overlapInd)

  # Window update
  phi_list <- update_windows(w_list, phiBar_list, model_list, 0, nT)

  # Consensus update
  phiBar_list <- update_consensus(phi_list, windows, overlapInd)

  # Dual update
  w_list <- mapply(update_dual, w_list, phi_list, phiBar_list,
                   MoreArgs = list(rho=rho), SIMPLIFY = FALSE)

  phiBar_k <- get_phiBar(phiBar_list, windows)

  dual_norm <- double(max_iter)
  primal_norm <- double(max_iter)

  iter <- 1

  while(iter <= max_iter){

    # Window update
    phi_list <- update_windows(w_list, phiBar_list, model_list, rho, nT)
    # Consensus update
    phiBar_list <- update_consensus(phi_list, windows, overlapInd)

    # Dual update
    w_list <- mapply(update_dual, w_list, phi_list, phiBar_list,
                     MoreArgs = list(rho=rho), SIMPLIFY = FALSE)
    
    # Convergence Metrics
    phiBar <- get_phiBar(phiBar_list, windows) 
    dual_resid <- phiBar - phiBar_k
    dual_norm[iter] <- 2*rho*Matrix::norm(dual_resid, type = "F")/(n_windows*overlap)
    phiBar_k <- phiBar

    primal_resid <- 0
    for (i in 1:n_windows){
      primal_resid <- primal_resid +
        norm(phi_list[[i]] - phiBar_list[[i]], "F")
    }

    primal_norm[iter] <- primal_resid/(n_windows*overlap)

    if (iter %% update == 0){
      print(sprintf("Iteration: %d Primal Resid Norm: %.4f Dual Resid Norm: %.4f", iter,
                    primal_norm[iter],
                    dual_norm[iter]))

    }

    if(dual_norm[iter] < eps*nT & primal_norm[iter] < eps*nT){
      primal_norm <- primal_norm[1:iter]
      dual_norm <- dual_norm[1:iter]
      print(sprintf("Converged in %d iterations", iter))
      break
    }
    iter <- iter+1
  }
  
  theta <- y - phiBar
  return(list(theta = theta, phiBar = phiBar, phi = phi_list,
              primal_norm = primal_norm, dual_norm = dual_norm,
              iter = iter))
}

