
gurobi_lp <- function(y, tau, lambda, D){
  n <- length(y)
  m <- nrow(D)
  
  model <- list()
  model$obj <- c(rep(tau, n), rep((1-tau), n), rep(lambda, 2*m))
  
  # Constraint Matrix
  model$A  <- as.matrix(cbind(D, -D, diag(m), -diag(m)))
  model$rhs <- as.numeric(D%*%y)
  # Constraint Types
  model$sense <- '='
  model$modelsense <- "min"
  
  # Type of variables (binary)
  params <- list(OutputFlag=0)
  result <- gurobi(model, params)
  print('Solution:')
  print(result$objval)
  
  theta <- y - result$x[1:n] + result$x[(n+1):(2*n)]
  
  return(theta)
}
