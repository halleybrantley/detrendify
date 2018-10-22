
gurobi_lp <- function(y, tau, lambda, D){
  n <- length(y)
  m <- nrow(D)
  np <- 2*n + 2*m
  nT <- length(tau)
  model <- list()
  
  for (i in 1:nT){
    model$obj <- c(model$obj, 
                   rep(tau[i], n), rep((1-tau[i]), n), rep(lambda, 2*m))
  
    model$rhs <- c(model$rhs, as.numeric(D%*%y))
  }
  
  model$rhs <- c(model$rhs, rep(0, nT-1))
  
  # Constraint Matrix
  if (length(tau) == 1){
    model$A  <- cbind(D, -D, diag(m), -diag(m))
  } else {
    model$A <- matrix(0, nrow =  (m+1)*nT - 1, ncol= np*nT ) 
    for (i in 1:nT){
      if (i < nT){
        model$A[m*nT+i, (1+(i-1)*np):(np*(i+1))] <- 
          c(rep(1, n), rep(-1, n), rep(0, 2*m),
            rep(-1, n), rep(1, n), rep(0, 2*m))
      }
      
      model$A[(1+m*(i-1)):(m*i), (1+np*(i-1)):(np*i)] <- 
                       cbind(D, -D, diag(m), -diag(m))
      
    }
  }
  
  # Constraint Types
  model$sense <- c(rep('=', m*nT), rep(">", nT-1))
  model$modelsense <- "min"
  
  # Type of variables (binary)
  params <- list(OutputFlag=0)
  result <- gurobi(model, params)
  print('Solution:')
  print(result$objval)
  
  theta <- matrix(0, nrow=n, ncol = nT)
  for (i in 1:nT){
    theta[,i] <- y - result$x[(1+np*(i-1)):(n+np*(i-1))] + 
      result$x[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }
  return(theta)
}
