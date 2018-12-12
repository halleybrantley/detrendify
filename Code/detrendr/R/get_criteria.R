checkloss <- function(e, tau){
  if (ncol(e) != length(tau)){
    stop("Number of columns in y must be same as length of tau")
  }
  obj <- e
  for (i in 1:length(tau)){
    obj[,i] <- obj[,i]*tau[i]
    obj[e[,i] < 0,i] <- e[e[,i] < 0,i]*(tau[i]-1)
  }
  return(obj)
  
}

get_criteria <- function(criteria, f_trend, y, tau, 
                         D = NULL, df_tol=1e-9, gamma=1, 
                         validID = NULL, yValid = NULL){
  
  if (criteria == "eBIC" || criteria == "SIC") {
    n <- length(y)
    resid_trend <- checkloss(y-f_trend, tau)
    df <- Matrix::colSums(abs(D%*%f_trend) > df_tol) 
    if (criteria == "eBIC"){
      scale_param <- 0.5 - abs(0.5-tau)
      BIC <- 2*colSums(resid_trend)/scale_param + log(n)*df +
        2*gamma*log(choose(nrow(D), df))
    } else {
      BIC <- log(colMeans(resid_trend)) + log(n)*df/(2*n)
    }
  } else if (criteria == "valid"){
    BIC <- colMeans(checkloss(yValid-f_trend[validID,,drop=FALSE], tau))
    df <- NA
  } 
  
  return(list(BIC = BIC, df=df))
} 