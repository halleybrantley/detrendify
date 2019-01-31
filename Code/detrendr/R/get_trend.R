#' Get quantile trends
#'
#' \code{get_trend} Returns the quantile trend matrix
#'
#' @param y observed data, should be equally spaced, may contain NA
#' @param tau quantile levels at which to evaluate trend
#' @param lambda penalty paramter controlling smoothness
#' @param k order of differencing
#' @importFrom utils installed.packages
#' @export
get_trend <- function(y, tau, lambda, k){
  mean_y <- mean(y, na.rm=T)
  sd_y <- sd(y, na.rm=T)
  y <- as.numeric(scale(y))*200
  model <- get_model(y, tau, lambda, k)
  pkgs <- installed.packages()[,"Package"]
  if("gurobi" %in% pkgs){
    solver <- "gurobi"
  } else if ("Rglpk" %in% pkgs){
    solver <- "Rglpk"
  } else {
    solver <- "lpSove"
  }
  theta <- solve_model(model, solver, y)
  return(theta/200*sd_y + mean_y)
}

