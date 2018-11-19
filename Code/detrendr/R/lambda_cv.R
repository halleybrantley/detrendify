#' Choose smoothing parameter using crossvalidation
#'
#' \code{lambda_cv} Selects smoothing parameter using cross validation
#'
#' @param y observed data
#' @param k order of differencing
#' @param tau quantile levels at which to evaluate trend
#' @param lambdaSeq smoothing penalty parameter options to compare
#' @param numFolds number of cross-validation folds to use
#' @param parallel TRUE/FALSE run in lambdas in parallel
#' @export
lambda_cv_kfold <- function(y, tau, k,
                      lambdaSeq = c(length(y)/10, length(y), length(y)*10),
                      numFolds = 5,
                      parallel = TRUE, 
                      cl = NULL){
  
  tau <- sort(tau)
  folds <- sample(1:numFolds, length(y), replace = TRUE)

  if (parallel){
    loss <- parLapply(cl, lambdaSeq, valid_checkLoss, folds,
                      y=y, tau=tau, k=k)
  } else {
    loss <- lapply(lambdaSeq, valid_checkLoss, folds,
                   y=y, tau=tau, k=k)
  }
  
  loss <- matrix(unlist(loss), ncol =length(tau), byrow = T)
  total <- rowSums(loss)
  return(list(lambda = lambdaSeq[which.min(total)],
              loss = loss))
}


#' Choose smoothing parameter leaving out every qth
#'
#' \code{lambda_valid} Selects smoothing parameter leaving out every qth observation
#'
#' @param y observed data
#' @param tau quantile levels at which to evaluate trend
#' @param k order of differencing
#' @param q observation to leave out
#' @param lambdaSeq smoothing penalty parameter options to compare
#' @export
lambda_valid <- function(y, tau, k, q,
                       lambdaSeq = seq(length(y)/10, length(y)*10, length(y)/10), 
                       df_tol = 1e-9){
  
  validID <- seq(q, length(y), q)
  yValid <- y[validID]
  y[validID] <- NA
  
  valid_checkloss <- c()
  
  for (i in 1:length(lambdaSeq)){
    lam <- lambdaSeq[i]
    f_trend <- gurobi_trend(y, tau, lam, k)
    valid_checkloss[i] <- mean(checkloss(yValid-f_trend[validID], tau))
  }
  
  plot(valid_checkloss~lambdaSeq, type="l", col="red")
  abline(v=lambdaSeq[which.min(valid_checkloss)], col="red")
  abline(h=min(valid_checkloss), col="red")
  lam <- lambdaSeq[which.min(valid_checkloss)]
  
  return(list(lambda = lam, 
              valid_checkloss = valid_checkloss))
}