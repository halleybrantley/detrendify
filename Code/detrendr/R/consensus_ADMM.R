
#' Consensus ADMM for overlapping windows
#'
#' \code{consensus ADMM} Uses ADMM to smooth overlapping windows
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param lambda smoothing penalty parameter
#' @param k order of differencing
#' @param rho step size for ADMM
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
#' w1 <- rep(0, overlap)
#' w2 <- rep(0, overlap)
#' z <- rep(0, overlap)
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
consensus_ADMM <- function(y, tau, lambda, k, rho, overlap, max_iter){
  y_n <- length(y)
  window_size <- floor((y_n+overlap)/2)

  Ind1 <- (window_size-overlap+1):(window_size)
  Ind2 <- 1:overlap
  
  y1 <- y[1:window_size]
  y2 <- y[(y_n-window_size+1):y_n]
  
  # Initial values
  w1 <- rep(0, length(y1))
  w2 <- rep(0, length(y2))
  z1 <- rep(0, length(y1))
  z2 <- rep(0, length(y2))
  rho <- 1

  
  # x-update
  phi1 <- quad_update(y1, tau, lambda, k, overlap, w1, z1, 0, first = TRUE)
  phi2 <- quad_update(y2, tau, lambda, k, overlap, w2, z2, 0, first = FALSE)
  
  # Average update
  phiBar <- numeric(length(y))
  phiBar[1:window_size] <- phi1 
  phiBar[(y_n-window_size+1):y_n] <- phi2
  phiBar[(y_n-window_size+1):(y_n-window_size+overlap)] <- 
    (phi1[Ind1] + phi2[Ind2])/2
  
  # Dual update
  w1 <- w1 + rho*(phi1 - phiBar[1:window_size])
  w2 <- w2 + rho*(phi2 - phiBar[(y_n-window_size+1):y_n])
  phiBar_k <- phiBar
  
  dual_norm <- double(max_iter)
  primal_norm <- double(max_iter)
  iter <- 1

  while(iter <= max_iter){
    phi1 <- quad_update(y1, tau, lambda, k, overlap, w1, z, rho, first = TRUE)
    phi2 <- quad_update(y2, tau, lambda, k, overlap, w2, z, rho, first = FALSE)
    
    # Average update
    phiBar <- numeric(length(y))
    phiBar[1:window_size] <- phi1 
    phiBar[(y_n-window_size+1):y_n] <- phi2
    phiBar[(y_n-window_size+1):(y_n-window_size+overlap)] <- 
      (phi1[Ind1] + phi2[Ind2])/2
    
    # Dual update
    w1 <- w1 + rho*(phi1 - phiBar[1:window_size])
    w2 <- w2 + rho*(phi2 - phiBar[(y_n-window_size+1):y_n])
    
    # Convergence Metrics
    dual_resid <- as.matrix(phiBar - phiBar_k)
    dual_norm[iter] <- 2*rho*Matrix::norm(dual_resid, type = "F")
    phiBar_k <- phiBar
    
    primal_resid <- 
      as.matrix(c(phi1 - phiBar[1:window_size],
                  phi2 - phiBar[(y_n-window_size+1):y_n]))
    
    primal_norm[iter] <- Matrix::norm(primal_resid, type = "F") 
    
    if(dual_norm[iter] < 0.0005 & primal_norm[iter] < 0.05){
      primal_norm <- primal_norm[1:iter]
      dual_norm <- dual_norm[1:iter]
      break
    }
    iter <- iter+1
  }
  
  theta <- y - phiBar
  return(list(theta = theta, phiBar = phiBar, phi1 = phi1, phi2 = phi2, 
              primal_norm = primal_norm, dual_norm = dual_norm))
}


#' Solve quandratic program step in consensus ADMM
#'
#' \code{quand_update} Uses gurobi to solve quadratic program step
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param lambda smoothing penalty parameter
#' @param k order of differencing
#' @param rho step size for ADMM
#' @param overlap integer length of overlap between windows
#' @param w value of dual variable
#' @param rho step size for ADMM
#' @param first boolean indicator of whether the first or second of the 
#' overlapping windows
quad_update <- function(y, tau, lambda, k, w, z, rho, first){
  tau <- sort(tau)
  D <- as.matrix(get_Dk(length(y), k))
  n <- length(y)
  m <- nrow(D)
  np <- 2*n + 2*m
  nT <- length(tau)
  missInd <- which(is.na(y))
  y[missInd] <- 0
  model <- list()
  
  for (i in 1:nT){
    model$obj <- c(model$obj,
                   rep(tau[i], n), rep((1-tau[i]), n), rep(lambda, 2*m))
    
    overlapInd <- (1+np*(i-1)):(2*n + np*(i-1))
    
    # Add in dual/consensus variables
    model$obj[overlapInd] <- model$obj[overlapInd] + c((w-rho*z), -(w-rho*z))

    model$obj[missInd] <- 0
    model$obj[missInd + n] <- 0
    model$rhs <- c(model$rhs, as.numeric(D%*%y))
  }
  
  model$rhs <- c(model$rhs, rep(0, n*(nT-1)))
  model$Q <- matrix(0, length(model$obj), length(model$obj))
  
  
  # Constraint Matrix
  if (length(tau) == 1){
    model$A  <- cbind(D, -D, diag(m), -diag(m))
    diag(model$Q[overlapInd, overlapInd]) <- rho/2
    
    
  } else {
    model$A <- matrix(0, nrow =  m*nT + n*(nT-1), ncol= np*nT )
    for (i in 1:nT){
      
      # Quadratic Term
      overlapInd <- (1+np*(i-1)):(2*n + np*(i-1))
      diag(model$Q[overlapInd, overlapInd]) <- rho/2
      
      # D%*%theta = eta constraint
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
  }
  
  # Constraint Types
  model$sense <- c(rep('=', m*nT), rep(">", n*(nT-1)))
  model$modelsense <- "min"
  
  # Type of variables (binary)
  params <- list(OutputFlag=0)
  result <- gurobi(model, params)
  
  phi <- matrix(0, nrow=n, ncol = nT)
  for (i in 1:nT){
    phi[,i] <- result$x[(1+np*(i-1)):(n+np*(i-1))] -
      result$x[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }
  
  return(phi)
}

