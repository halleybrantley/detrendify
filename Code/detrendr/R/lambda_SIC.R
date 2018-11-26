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
                       lambdaSeq = seq(length(y)/5, length(y)*5, length(y)/5), 
                       df_tol = 1e-9, 
                       plot_lambda = FALSE){
  
  df_trend <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
  SIC_trend <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
  
  for (i in 1:length(lambdaSeq)){
    lam <- lambdaSeq[i]
    f_trend <- gurobi_trend(y, tau, lam, k)
    resid_trend <- checkloss(y-f_trend, tau)
    df_trend[i,] <- colSums(resid_trend<df_tol)
    SIC_trend[i,] <- log(colMeans(resid_trend)) + log(n)*df_trend[i,]/(2*n)
  }
  
  SIC_scale <- as.data.frame(scale(SIC_trend))
  SIC_scale$mean <- rowSums(SIC_scale)
  
  if (plot_lambda){
    plot(SIC_scale[,1]~lambdaSeq, type="l", col="red", 
         ylim = c(min(SIC_scale), max(SIC_scale)))
    for (i in 1:length(tau)){
      lines(SIC_scale[,i]~lambdaSeq)
    }
    lines(SIC_scale$mean ~ lambdaSeq, col="red")
    abline(v=lambdaSeq[which.min(SIC_scale$mean)], col="red")
    abline(h=min(SIC_scale$mean), col="red")
  }

  lam <- lambdaSeq[which.min(SIC_scale$mean)]

  
  return(list(lambda = lam, 
              SIC = SIC_trend, 
              df = df_trend))
}