#' Get quantile trends
#'
#' \code{get_trend} Returns the quantile trend matrix
#'
#' @param y observed data, should be equally spaced, may contain NA
#' @param tau quantile levels at which to evaluate trend
#' @param lambda penalty paramter controlling smoothness
#' @param k order of differencing
#' @export
get_trend <- function(y, tau, lambda, k){
  model <- get_model_data(y, tau, lambda, k)
  pkgs <- installed.packages()[,"Package"]
  if("gurobi" %in% pkgs){
    solver <- "gurobi"
  } else if ("Rglpk" %in% pkgs){
    solver <- "Rglpk"
  } else {
    solver <- "lpSove"
  }
  theta <- solve_model(model, y, solver)
  return(theta)
}

solve_model <- function(model, y, solver){
  if (solver == "gurobi"){
    require(gurobi)
    theta <- solve_gurobi(model, y)
  } else if (solver == "Rglpk") {
    require(Rglpk)
    theta <- solve_glpk(model, y)
  } else if (solver == "lpSolve"){
    theta <- solve_lp(model, y)
  } else {
    stop("Solver must be one of 'gurobi', 'Rglpk', or 'lpSolve'")
  }
  return(theta)
}  