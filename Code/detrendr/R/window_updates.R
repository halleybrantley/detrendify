#' Update windows
#'
#' \code{update_windows} Minimizes the Lagrangian for each window
#'
#' @param y_list
#' @param w_list
#' @param phiBar_list
#' @param tau
#' @param lambda
#' @param k
#' @param rho step size for ADMM
#' @export
update_windows <- function(w_list, phiBar_list, model_list, rho, nT){
  model_list <- mapply(update_model, model_list, w_list, phiBar_list, 
                       rho=rho, nT=nT, SIMPLIFY = FALSE)
  phi_list <- lapply(model_list, solve_model, solver="gurobi", trend = FALSE)
  return(phi_list)
}

#' Update consensus variable (overlap)
#'
#' \code{update_consensus}
#'
#' @param phi_list
#' @param windows
#' @param overlapInd
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
#' @param phi_list
#' @param windows
#' @param overlapInd
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
#' @param w
#' @param phi
#' @param phiBar
#' @param rho
#' @export
update_dual <- function(w, phi, phiBar, rho) {
  w + as.numeric(rho*(phi - phiBar))
}

#' Update model with ADMM step parameter and dual and consensus variables
#'
#' \code{update_dual}
#'
#' @param w
#' @param phi
#' @param phiBar
#' @param rho
#' @export
update_model <- function(model, w, z, rho, nT){
  n <- length(w)/nT
  np <- length(model$obj)/nT
  
  for (i in 1:nT){
    thetaInd <- (1+np*(i-1)):(2*n + np*(i-1))
    dualInd <- (1 + n*(i-1)):(n + n*(i-1))
    # Add in dual/consensus variables
    model$obj[thetaInd] <- model$obj[thetaInd]  +
      c((w[dualInd]-rho*z[dualInd]), -(w[dualInd]-rho*z[dualInd]))
  }
  if (is.null(model$Q)) {
    thetaAll <- c()
    for (i in 1:nT){
      thetaInd <- (1+np*(i-1)):(2*n + np*(i-1))
      thetaAll <- c(thetaAll, thetaInd)
    }
    model$Q <- sparseMatrix(thetaAll, thetaAll, x=rho/2, 
                            dims = c(length(model$obj), length(model$obj)))
  }
  return(model)
}


