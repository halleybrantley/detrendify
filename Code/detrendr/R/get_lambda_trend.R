#' Choose smoothing parameter using SIC
#'
#' \code{lambda_SIC} Selects smoothing parameter using Schwarz Information
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param k order of differencing
#' @param lambdaSeq smoothing penalty parameter options to compare
#' @export
get_lambda_trend <- function(y, tau, k,
                       lambdaSeq = length(y)^seq(0, 1.5, length.out=20), 
                       df_tol = 1e-9, 
                       gamma = 1,
                       plot_lambda = FALSE, 
                       solver = NULL, 
                       criteria = "eBIC"){
  
  if(!(criteria %in% c("eBIC", "valid", "SIC"))){
    stop("criteria must be one of 'eBIC', 'valid', 'SIC'")
  }
  
  # Set linear program solver
  if (is.null(solver)){
    pkgs <- installed.packages()[,"Package"]
    if("gurobi" %in% pkgs){
      solver <- "gurobi"
    } else if ("Rglpk" %in% pkgs){
      solver <- "Rglpk"
    } else {
      solver <- "lpSove"
    }
  }
  
  if (criteria == "valid"){
    validID <- seq(5, length(y), 5)
    yValid <- y[validID]
    y[validID] <- NA
    D <- NULL
  } else {
    validID <- NULL
    yValid <- NULL
    D <- get_Dk(n, k)
  }
  
  df <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
  BIC <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))

  model <- get_model_data(y, tau, lambdaSeq[1], k)
  n <- length(y)
  m <- n-k
  missInd <- which(is.na(y))
  
  for (i in 1:length(lambdaSeq)){
    model$obj <- get_obj(tau, rep(lambdaSeq[i], length(tau)), n, m, missInd)
    f_trend <- solve_model(model, y, solver)
    model_crit <- get_criteria(criteria, f_trend, y, tau, 
                               D, df_tol, gamma, 
                               validID, yValid)
    df[i,] <- model_crit$df
    BIC[i,] <- model_crit$BIC
  }
  
  lambda <- lambdaSeq[apply(BIC, 2, which.min)]
  
  if (plot_lambda){
    plot(BIC[,1]~log(lambdaSeq), type="l", col="red", 
         ylim = c(min(BIC), max(BIC)))
    for (i in 1:length(tau)){
      lines(BIC[,i]~log(lambdaSeq))
    }
    abline(v=log(lambda))
  }
  
  model$obj <- get_obj(tau, lambda, n, m, missInd)
  theta <- solve_model(model, y, solver)

  return(list(theta = theta,
              lambda = lambda, 
              BIC = BIC, 
              df = df))
}