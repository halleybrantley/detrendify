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
#' @param eps 
#' @param update number of iterations at which to print residuals
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
get_trend_windows <- function(y, tau, lambda, k, rho=1, window_size,
                           overlap, max_iter, update=10, 
                           quad = TRUE, use_gurobi = TRUE, 
                           eps_abs = 1e-2, 
                           eps_rel = 1e-3){
  if (use_gurobi){
    solver <- "gurobi"
  } else {
    # First estimate uses LP not QP
    solver <- "lpSolve"
  }
  
  window_size <- round(window_size)
  
  y_n <- length(y)
  tau <- sort(tau)
  nT <- length(tau)
  D <- get_Dk(window_size, k)
  n_windows <- ceiling((y_n-overlap)/(window_size-overlap))
  windows <- matrix(FALSE, y_n, n_windows)
  y_list <- list()
  w_list <- list()
  
  # Initial values
  for (i in 1:n_windows){
    start_I <- 1+(window_size-overlap)*(i-1)
    end_I <- min((window_size + (window_size-overlap)*(i-1)), y_n)
    windows[start_I:end_I,i] <- TRUE
    y_list[[i]] <- y[start_I:end_I]
    len <- end_I - start_I + 1
    w_list[[i]] <- rep(0, (2*len - k)*length(tau))
  }
  overlapInd <- rowSums(windows) > 1
  
  # Window initial LP fit
  model_list <- lapply(y_list, get_model, tau=tau, lambda=lambda, k=k)
  phi_list <- mapply(solve_model, model_list, y_list, 
                     solver = solver, trend=FALSE, SIMPLIFY = FALSE)
  eta_list <- mapply(get_eta, phi_list, y_list, k=3, SIMPLIFY = FALSE)
  
  # Change to QP solver
  if (solver == "lpSolve"){
    solver <- "quadprog"
  }
  # Consensus update
  phiBar_list <- update_consensus(phi_list, windows, overlapInd)
  etaBar_list <- mapply(get_eta, phiBar_list, y_list, k=k, SIMPLIFY = FALSE)
  
  num_param <- sum(sapply(phiBar_list, length))
  
  if (quad){
    # Use w_list for z since it is all zeros, don't want to store current z
    eta0 <- etaBar_list
    phi0 <- phiBar_list
    for (i in 1:length(eta0)){
      eta0[[i]][] <- 0
      phi0[[i]][] <- 0
    }
    model_list <- mapply(update_model, model_list, w_list, phi0,
                         eta0,
                         rho=rho, nT=nT, SIMPLIFY = FALSE)
  } else {
    model_list <- mapply(get_model_abs, model_list, phiBar_list, rho, nT, 
                         SIMPLIFY = FALSE)  
  }
  
  # Dual update
  w_list <- mapply(update_dual, w_list, 
                   phi_list, phiBar_list, 
                   eta_list, etaBar_list, 
                   MoreArgs = list(rho=rho), SIMPLIFY = FALSE)

  dual_norm <- double(max_iter)
  primal_norm <- double(max_iter)
  phiBar_listk <- phiBar_list
  iter <- 1

  while(iter <= max_iter){

    # Window update
    phi_list <- update_windows(w_list, phiBar_list, etaBar_list, 
                               model_list, rho, nT, 
                               quad, solver)
    # Consensus update
    phiBar_list <- update_consensus(phi_list, windows, overlapInd)
    etaBar_list <- mapply(get_eta, phiBar_list, y_list, k=k, SIMPLIFY = FALSE)
    
    # Dual update
    w_list <- mapply(update_dual, w_list, 
                     phi_list, phiBar_list, 
                     eta_list, etaBar_list, 
                     MoreArgs = list(rho=rho), SIMPLIFY = FALSE)
    
    # Convergence Metrics

    dual_norm[iter] <- sqrt(rho*sum(mapply(list_diff_norm, phiBar_list, phiBar_listk)))
    primal_norm[iter] <- sqrt(sum(mapply(list_diff_norm, phi_list, phiBar_list)))
    eps_pri <- sqrt(num_param)*eps_abs + 
      eps_rel*max(sqrt(c(sapply(phi_list, Matrix::norm, type="F"), 
                    sapply(phiBar_list, Matrix::norm, type="F"))))
    
    eps_dual <- sqrt(num_param)*eps_abs + 
      eps_rel*sqrt(sum(sapply(w_list, norm, type="2")))
    
    phiBar_listk <- phiBar_list
    
    if (iter %% update == 0){
      print(sprintf("Iteration: %d Primal Resid Norm: %.4f eps_pri: %.4f, Dual Resid Norm: %.4f  eps_dual %.4f", 
                    iter,
                    primal_norm[iter],
                    eps_pri,
                    dual_norm[iter],
                    eps_dual))

    }

    if (iter == 1 ){
      primal_resid0 <-  primal_norm[iter] 
      dual_resid0 <- dual_norm[iter]
    } else if(dual_norm[iter] < eps_dual & 
              primal_norm[iter] < eps_pri){
      primal_norm <- primal_norm[1:iter]
      dual_norm <- dual_norm[1:iter]
      print(sprintf("Converged in %d iterations", iter))
      break
    }
    iter <- iter+1
  }
  y[is.na(y)] <- 0
  theta <- y - get_phiBar(phiBar_list, windows)
  return(theta)
}

list_diff_norm <- function(list1, list2) {
  return(Matrix::norm(list1 - list2, type = "F"))
}
