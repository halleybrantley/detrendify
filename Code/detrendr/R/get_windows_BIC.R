#' Choose smoothing parameter using SIC
#'
#' \code{get_windows_BIC} Selects smoothing parameter using Schwarz Information
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param k order of differencing
#' @param lambdaSeq smoothing penalty parameter options to compare
#' @param df_tol tolerance for determining degrees of freedom (D%*%theta > df_tol)
#' @param parameter for eBIC
#' @param plot_lambda TRUE/FALSE for plotting lambda by model criteria
#' @param solver LP solver, can be "gurobi", "Rglpk", or "lpSolve"
#' @param criter criteria to use for lambda selection, must be "eBIC", "SIC", or 
#' "valid"  
#' @export
get_windows_BIC <- function(y, tau, k, window_size, overlap,
                             lambdaSeq = length(y)^seq(0, 1.5, length.out=20), 
                             df_tol = 1e-9, 
                             gamma = 1,
                             plot_lambda = FALSE, 
                             solver = NULL, 
                             criteria = "eBIC", 
                             max_iter = 10){
  
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

  n <- length(y)
  m <- n-k
  missInd <- which(is.na(y))  
  
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

  
  for (i in 1:length(lambdaSeq)){
    f_trend <- get_trend_windows(y, tau, lambdaSeq[i], k=k, window_size = window_size,
                                 overlap=overlap, max_iter=max_iter, update=5, 
                                 quad = TRUE)
    model_crit <- get_criteria(criteria, f_trend, y, tau, 
                               D, df_tol, gamma, 
                               validID, yValid)
    df[i,] <- model_crit$df
    BIC[i,] <- model_crit$BIC
    print(sprintf("i=%d lambda=%f", i, lambdaSeq[i]))
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
  
  f_trend <- get_trend_windows(y, tau, lambda, k, window_size,
                               overlap, max_iter=max_iter, update=5, 
                               quad = TRUE)
  
  return(list(trend = f_trend,
              lambda = lambda, 
              BIC = BIC, 
              df = df))
}
