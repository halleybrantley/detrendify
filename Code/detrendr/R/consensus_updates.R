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
  n_windows <- length(w_list)
  phi_list <- list()
  #foreach(i = 1:n_windows) %dopar% {
  for (i in 1:n_windows){
  phi_list[[i]] <- quad_update(w_list[[i]], as.numeric(phiBar_list[[i]]), 
                model_list[[i]], rho, nT)
  }
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



create_model <- function(y, tau, lambda, D){
  n <- length(y)
  missInd <- which(is.na(y))
  y[missInd] <- 0
  
  if(ncol(D) != length(y)){
    D <- D[1:(length(y) - ncol(D)+nrow(D)), 1:length(y)]
  }
  m <- nrow(D) 
  model <- list()
  np <- 2*n + 2*m
  nT <- length(tau)

  for (i in 1:nT){
    model$obj <- c(model$obj,
                   rep(tau[i], n), rep((1-tau[i]), n), rep(lambda, 2*m))
    model$rhs <- c(model$rhs, as.numeric(D%*%y))
  }
  model$obj[missInd] <- 0
  model$obj[missInd + n] <- 0
  
  model$rhs <- c(model$rhs, rep(0, n*(nT-1)))
  model$Q <- matrix(0, length(model$obj), length(model$obj))
  
  # Constraint Matrix
  model$A <- matrix(0, nrow =  m*nT + n*(nT-1), ncol= np*nT )
  for (i in 1:nT){
    # Quadratic Term
    thetaInd <- (1+np*(i-1)):(2*n + np*(i-1))
    diag(model$Q[thetaInd, thetaInd]) <- rho/2
    
    # Constraint for D%*%theta = eta 
    model$A[(1+m*(i-1)):(m*i), (1+np*(i-1)):(np*i)] <-
      cbind(D, -D, diag(m), -diag(m))
    
    # Non-crossing theta(tau) constrains
    if (i < nT){
      model$A[(m*nT+1+n*(i-1)):(m*nT+n*(i)), (1+(i-1)*np):(2*n + (i-1)*np)] <-
        cbind(diag(n), -diag(n))
      
      model$A[(m*nT+1+n*(i-1)):(m*nT+n*(i)),
              (1+i*np):(2*n + i*np)] <-
        cbind(-diag(n), diag(n))
    }
  }
  
  # Constraint Types
  model$sense <- c(rep('=', m*nT), rep(">", n*(nT-1)))
  model$modelsense <- "min"
  
  return(model)
  
}


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
  if (model$Q[1,1] != rho/2) {
    model$Q[model$Q != 0] <- rho/2
  }
  return(model)
}

#' Solve quadratic program step in consensus ADMM
#'
#' \code{quad_update} Uses gurobi to solve quadratic program step
#'
#' @param w value of dual variable
#' @param z consensus variable
#' @param rho step size for ADMM
#' @param nT number of quantiles
#' @export
quad_update <- function(w, z, model, rho, nT){
  if(length(w) != length(z)) {
    stop("w and z must be the same length")
  }

  params <- list(OutputFlag=0)
  model <- update_model(model, w, z, rho, nT)
  result <- gurobi(model, params)
  np <- length(model$obj)/nT
  n <- length(w)/nT
  
  phi <- matrix(0, nrow=n, ncol = nT)
  for (i in 1:nT){
    phi[,i] <- result$x[(1+np*(i-1)):(n+np*(i-1))] -
      result$x[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }

  return(phi)
}
