#' Solves the model using chosen solver
#'
#' \code{solve_model} Returns the trend
#'
#' @param y observed data, should be equally spaced, may contain NA
#' @param tau quantile levels at which to evaluate trend
#' @param lambda penalty paramter controlling smoothness
#' @param k order of differencing
#' @param trend if TRUE returns trend, if FALSE returns residuals
#' @export

solve_model <- function(model, solver, y=NULL, trend = TRUE){
  if (trend && is.null(y)){
    stop("y required to get trend")
  }
  
  if (solver == "gurobi"){
    require(gurobi)
    x <- solve_gurobi(model)
  } else if (solver == "Rglpk") {
    require(Rglpk)
    x <- solve_glpk(model)
  } else if (solver == "lpSolve"){
    x <- solve_lp(model)
  } else {
    stop("Solver must be one of 'gurobi', 'Rglpk', or 'lpSolve'")
  }
  
  nT <- model$nT
  n <- model$n
  np <- model$np
  phi <- matrix(0, nrow=n, ncol = nT)
  
  for (i in 1:nT){
    phi[,i] <- x[(1+np*(i-1)):(n+np*(i-1))] - 
      x[(n+1 + np*(i-1)):(2*n + np*(i-1))]
  }

  if (trend){
    y[which(is.na(y))] <- 0
    return(y-phi)
  } else {
    return(phi)
  }
}  


# Wrappers for different solvers

solve_gurobi <- function(model){
  params <- list(OutputFlag=0)
  result <- gurobi(model, params)
  iter <- 1
  while(result$status != "OPTIMAL"){
    print("Problem not solved, adding 0.01 to lambda.")
    lambda <- model$obj[length(model$obj)]
    model$obj[model$obj == lambda] <- lambda + 0.01
    iter <- iter+1
    if(iter == 10){
      stop("Problem not solved.")
    }
  }
  return(result$x)
}

solve_glpk <- function(model_list){
  result <- Rglpk_solve_LP(obj = model_list$obj, 
                           mat = model_list$A, 
                           dir = paste0(model_list$sense, "="), 
                           rhs = model_list$rhs)
  return(result$solution)
}

solve_lp <- function(model_list){
  result <- lp(objective.in = model_list$obj, 
               const.mat = as.matrix(model_list$A), 
               const.dir = paste0(model_list$sense, "="), 
               const.rhs = model_list$rhs)
  return(result$solution)
}

