#' Choose smoothing parameter using SIC
#'
#' \code{lambda_SIC} Selects smoothing parameter using Schwarz Information
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param k order of differencing
#' @param lambdaSeq smoothing penalty parameter options to compare
#' @export
lambda_eBIC <- function(y, tau, k,
                       lambdaSeq = length(y)^seq(0, 2, length.out=50), 
                       df_tol = 1e-9, 
                       gamma = 1,
                       plot_lambda = FALSE, 
                       single_lambda = FALSE){
  
  df_trend <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
  SIC_trend <- matrix(NA, nrow=length(lambdaSeq), ncol=length(tau))
  tau_min <- sapply(tau, function(x){min(x,1-x)})
  D <- get_Dk(length(y), k)
  
  for (i in 1:length(lambdaSeq)){
    lam <- lambdaSeq[i]
    suppressMessages(f_trend <- gurobi_trend(y, tau, lam, k))
    resid_trend <- checkloss(y-f_trend, tau)
    discr_diff <- abs(D%*%f_trend)
    df_trend[i,] <- Matrix::colSums(discr_diff > df_tol) #colSums(abs(y-f_trend)<df_tol)
    SIC_trend[i,] <- 2*colSums(resid_trend) + tau_min*log(n)*df_trend[i,] +
      2*gamma*log(choose(n-k, df_trend[i,]))
  }
  
  SIC_scale <- as.data.frame(scale(SIC_trend))
  SIC_scale$mean <- rowMeans(SIC_scale)
  
  if (single_lambda){
    lam <- lambdaSeq[which.min(SIC_scale$mean)]
  } else {
    lam <- lambdaSeq[apply(SIC_trend, 2, which.min)]
  }  
  
  
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

  

  
  return(list(lambda = lam, 
              SIC = SIC_trend, 
              df = df_trend))
}