# Wrappers for different solvers

solve_gurobi <- function(model_list, y){
  params <- list(OutputFlag=0)
  result <- gurobi(model_list, params)
  nT <- model_list$nT
  n <- model_list$n
  np <- model_list$np
  theta <- matrix(0, nrow=n, ncol = nT)
  y[which(is.na(y))] <- 0
  for (i in 1:nT){
    theta[,i] <- y - result$x[(1+np*(i-1)):(n+np*(i-1))] +
      result$x[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }
  return(theta)
}

solve_glpk <- function(model_list, y){
  result <- Rglpk_solve_LP(obj = model_list$obj, 
                           mat = model_list$A, 
                           dir = paste0(model_list$sense, "="), 
                           rhs = model_list$rhs)
  nT <- model_list$nT
  n <- model_list$n
  np <- model_list$np
  theta <- matrix(0, nrow=n, ncol = nT)
  y[which(is.na(y))] <- 0
  for (i in 1:nT){
    theta[,i] <- y - result$solution[(1+np*(i-1)):(n+np*(i-1))] +
      result$solution[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }
  return(theta)
}

solve_lp <- function(model_list, y){
  result <- lp(objective.in = model_list$obj, 
               const.mat = as.matrix(model_list$A), 
               const.dir = paste0(model_list$sense, "="), 
               const.rhs = model_list$rhs)
  nT <- model_list$nT
  n <- model_list$n
  np <- model_list$np
  theta <- matrix(0, nrow=n, ncol = nT)
  y[which(is.na(y))] <- 0
  for (i in 1:nT){
    theta[,i] <- y - result$solution[(1+np*(i-1)):(n+np*(i-1))] +
      result$solution[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }
  return(theta)
}
