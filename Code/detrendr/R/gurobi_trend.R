#' Quantile Regression with Trend Filtering (Gurobi Solver)
#'
#' \code{gurobi_trend} Determines the quantile regression trend filtering solution
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param lambda penalty paramter controlling smoothness
#' @param k order of differencing
#' @export
gurobi_trend <- function(y, tau, lambda, k){

  if(tau >= 1 || tau <= 0){
    stop("tau must be between 0 and 1.")
  }

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
    model$obj[missInd] <- 0
    model$obj[missInd + n] <- 0
    model$rhs <- c(model$rhs, as.numeric(D%*%y))
  }

  model$rhs <- c(model$rhs, rep(0, n*(nT-1)))

  # Constraint Matrix
  if (length(tau) == 1){
    model$A  <- cbind(D, -D, diag(m), -diag(m))
  } else {
    model$A <- matrix(0, nrow =  m*nT + n*(nT-1), ncol= np*nT )
    for (i in 1:nT){

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

  theta <- matrix(0, nrow=n, ncol = nT)
  for (i in 1:nT){
    theta[,i] <- y - result$x[(1+np*(i-1)):(n+np*(i-1))] +
      result$x[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }
  return(theta)
}



