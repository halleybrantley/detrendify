#' Calculate total check loss for validation data
#'
#' \code{valid_checkLoss} Calculate checkloss for validation data
#'
#' @param lambda smoothing penalty parameter
#' @param folds vector indicating folds for k-fold cross-validation
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param k order of differencing
valid_checkLoss <- function(lambda, folds, y, tau, k){
  loss <- numeric(length(tau))
  
  for (fold in unique(folds)){
    train <- y
    holdout <- (folds == fold)
    train[holdout] <- NA
    theta_gurobi <- gurobi_trend(train, tau, lambda, k)
    
    for (j in 1:length(tau)){
      loss[j] <- loss[j] + check_loss(y[holdout]-theta_gurobi[holdout,j], tau[j])
    }
    
  }
  return(loss)
}
