# Functions for choosing smoothing parameter

#' Evaluate checkloss function
#'
#' \code{checkloss}
#'
#' @param e argument of checkloss function 
#' @param tau quantile to be used
#' @export
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

H <- function(x) {
  return(x*log(1/x) + (1-x)*log(1/(1-x)))
}

# Criteria for evaluating the smoothness parameter
#'
#' \code{get_criteria}
#'
#' @param criteria label of criteria to be used, must be one of "eBIC", "SIC", "valid"
#' @param f_trend fitted quantile trend(s)
#' @param observed data vector
#' @param tau quantile vector
#' @param D discrete differencing matrix (for SIC and eBIC)
#' @param df_tol tolerance for determining degrees of freedom (D%*%theta > df_tol)
#' @param gamma parameter for eBIC
#' @param validID index of data to be used for validation (for valid method) 
#' @param yValid validation data (for valid method) 
#' @export
get_criteria <- function(criteria, f_trend, y, tau, 
                         D = NULL, df_tol=1e-9, gamma=1, 
                         validID = NULL, yValid = NULL){
  
  if (criteria == "eBIC" || criteria == "SIC") {
    n <- length(y)
    missInd <- which(is.na(y))
    resid_trend <- checkloss(y[-missInd]-f_trend[-missInd,,drop=FALSE], tau)
    df <- Matrix::colSums(abs(D%*%f_trend) > df_tol) 
    if (criteria == "eBIC"){
      scale_param <- 0.5 - abs(0.5-tau)
      BIC <- 2*colSums(resid_trend, na.rm=T)/scale_param + log(n)*df +
        2*gamma*log(choose(nrow(D), df))
      if (is.infinite(BIC)) {
        BIC <- 2*colSums(resid_trend, na.rm=T)/scale_param + log(n)*df +
          2*gamma*nrow(D)*H(df/nrow(D))
      }
    } else {
      BIC <- log(colMeans(resid_trend, na.rm=T)) + log(n)*df/(2*n)
    }
  } else if (criteria == "valid"){
    missInd <- is.na(yValid)
    BIC <- colMeans(checkloss(
      yValid[-missInd]-f_trend[validID[-missInd],,drop=FALSE], tau))
    df <- NA
  } 
  
  return(list(BIC = BIC, df=df))
} 