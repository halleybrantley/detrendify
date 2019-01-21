#' Update windows
#'
#' \code{update_windows} Minimizes the Lagrangian for each window
#'
#' @param y_list Observed data
#' @param w_list Dual variable
#' @param phiBar_list Consensus variable
#' @param tau quantiles
#' @param lambda smoothing parameter
#' @param k order of differencing matrix
#' @param rho step size for ADMM
#' @export
update_windows <- function(w_list, phiBar_list, etaBar_list, 
                           model_list, rho, nT, 
                           quad=FALSE, solver="gurobi"){
  if (quad){
    model_list <- mapply(update_model, model_list, w_list, phiBar_list,
                         etaBar_list,
                         rho=rho, nT=nT, SIMPLIFY = FALSE)
  } else {
    model_list <- mapply(update_model_abs, model_list, w_list, phiBar_list, 
                         nT=nT, SIMPLIFY = FALSE)
  }
  phi_list <- lapply(model_list, solve_model, solver=solver, trend = FALSE)
  return(phi_list)
}

#' Update consensus variable (overlap)
#'
#' \code{update_consensus}
#'
#' @param phi_list Primal variables
#' @param windows Matrix indicating window group
#' @param overlapInd Indices of overlap between windows
#' @export
update_consensus <- function(phi_list, windows, overlapInd){

  phiBar <- matrix(0, nrow = nrow(windows), ncol = length(tau))
  n_windows <- length(phi_list)
  
  for (i in 1:n_windows){
    phiBar[windows[,i],] <- phiBar[windows[,i],] + phi_list[[i]]
  }
  phiBar[overlapInd, ] <- phiBar[overlapInd, ]/2

  phiBar_list <- list()
  for (i in 1:n_windows){
    phiBar_list[[i]] <- phiBar[windows[,i],]
  }
  return(phiBar_list)
}

#' Update consensus variable (overlap)
#'
#' \code{get_phiBar}
#'
#' @param phiBar_list List of consensus variable
#' @param windows Matrix indicating window group
#' @export
get_phiBar <- function(phiBar_list, windows){
  phiBar <- matrix(0, nrow = nrow(windows), ncol = length(tau))
  for (i in 1:length(phiBar_list)){
     phiBar[windows[,i],] <- phiBar_list[[i]] 
  }
  return(phiBar)
}


#' Update dual variable (overlap)
#'
#' \code{update_dual}
#'
#' @param w dual variable
#' @param phi primal variable
#' @param phiBar consensus variable
#' @param rho ADMM step size parameter
#' @export
update_dual <- function(w, phi, phiBar, eta, etaBar, rho) {
  w + as.numeric(rbind(rho*(phi - phiBar), rho*(eta-etaBar)))
}

#' Update model with ADMM step parameter and dual and consensus variables
#'
#' \code{update_model}
#'
#' @param model List of QP model elements
#' @param w dual variable
#' @param z Consensus variable
#' @param rho ADMM step size
#' @param nT Number of quantiles being estimated
#' @export
update_model <- function(model, w, phiBar, etaBar, rho, nT){
  z <- as.numeric(rbind(phiBar, etaBar))
  n <- length(phiBar)/nT
  np <- length(model$obj)/nT
  
  for (i in 1:nT){
    objInd <- (1+np*(i-1)):(np*i)
    dual_theta <- (1 + np/2*(i-1)):(np/2*(i-1) + n)
    dual_eta <- (1 + np/2*(i-1) + n):(np/2*i)
    # Add in dual/consensus variables
    model$obj[objInd] <- model$obj[objInd]  +
      c(w[dual_theta]-rho*z[dual_theta], -(w[dual_theta]-rho*z[dual_theta]),
        w[dual_eta]-rho*z[dual_eta], -(w[dual_eta]-rho*z[dual_eta]))
  }
  
  if (is.null(model$Q)) {
    model$Q <- Matrix::Diagonal(length(model$obj), rho/2)
  }
  return(model)
}


#' Update model with ADMM step parameter and dual and consensus variables
#'
#' \code{update_model_abs}
#'
#' @param model Model object
#' @param w dual variable
#' @param z consensus variable
#' @param nT number of quantile levels
#' @export
update_model_abs <- function(model, w, z, nT){
  n <- length(w)/nT
  np <- length(model$obj)/nT
  
  for (i in 1:nT){
    thetaInd <- (1+np*(i-1)):(2*n + np*(i-1))
    dualInd <- (1 + n*(i-1)):(n + n*(i-1))
    # Add in dual/consensus variables
    model$obj[thetaInd] <- model$obj[thetaInd]  +
      c(w[dualInd], -w[dualInd])
  }
  nC <- length(model$rhs)
  model$rhs[(nC-nT*n+1):nC] <- as.numeric(z)
  return(model)
}

#' Get model with ADMM step parameter and dual and consensus variables
#'
#' \code{get_model_abs}
#'
#' @param model Model object
#' @param z consensus variable
#' @param rho ADMM step parameter
#' @param nT number of quantile levels
#' @export
get_model_abs <- function(model, z, rho, nT){
  n <- length(z)/nT
  n0 <- length(model$obj)
  np <- n0/nT
  
  model$obj <- c(model$obj, rep(rho/2, 2*n*nT))
  newConst <- Matrix(0, nrow = n*nT, ncol= length(model$obj), sparse=TRUE)
  
  for (i in 1:nT){
    newConst[(1+n*(i-1)):(n*i),(n0+1+2*n*(i-1)):(n0+2*n*i)] <- 
      cbind(-Diagonal(n), Diagonal(n))
    newConst[(1+n*(i-1)):(n*i),(1+np*(i-1)):(2*n + np*(i-1))] <- 
      cbind(Diagonal(n), -Diagonal(n))
  }
  
  A <- Matrix(0, nrow = nrow(model$A), ncol = 2*n*nT, sparse = TRUE) 
  model$A <- rbind(cbind(model$A, A), newConst)
  model$rhs <- c(model$rhs, as.numeric(z))
  model$sense <- c(model$sense , rep("=", length(z)))
  return(model)
}


get_eta <- function(phi, y, k){
  D <- get_Dk(length(y), k)
  eta <- D%*%y - D%*%phi
  eta[is.na(eta)] <- 0
  return(as.matrix(eta))
}


