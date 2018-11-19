#' Choose smoothing parameter using SIC
#'
#' \code{lambda_SIC} Selects smoothing parameter using Schwarz Information
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param k order of differencing
#' @param lambdaSeq smoothing penalty parameter options to compare
#' @export
lambda_SIC <- function(y, tau, k,
                       lambdaSeq = seq(length(y)/10, length(y)*10, length(y)/10), 
                       df_tol = 1e-9){
  
  df_trend <- c()
  SIC_trend <- c()
  
  for (i in 1:length(lambdaSeq)){
    lam <- lambdaSeq[i]
    f_trend <- gurobi_trend(y, tau, lam, k)
    resid_trend <- checkloss(y-f_trend, tau)
    df_trend[i] <- sum(resid_trend<df_tol)
    SIC_trend[i] <- log(mean(resid_trend)) + log(n)*df_trend[i]/(2*n)
  }
  
  plot(SIC_trend~lambdaSeq, type="l", col="red")
  abline(v=lambdaSeq[which.min(SIC_trend)], col="red")
  abline(h=min(SIC_trend), col="red")
  lam <- lambdaSeq[which.min(SIC_trend)]

  
  return(list(lambda = lam, 
              SIC = SIC_trend, 
              df = df_trend))
}