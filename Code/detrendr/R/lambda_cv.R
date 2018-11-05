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
lambda_cv <- function(y, tau, k,
                      lambdaSeq = c(length(y)/10, length(y), length(y)*10),
                      numFolds = 5,
                      parallel = TRUE){
  
  tau <- sort(tau)
  folds <- sample(1:numFolds, length(y), replace = TRUE)
  no_cores <- detectCores() - 1
  
  if (parallel){
    cl <- makeCluster(no_cores, type="FORK")
    loss <- parLapply(cl, lambdaSeq, valid_checkLoss, folds,
                      y=y, tau=tau, k=k)
    stopCluster(cl)
  } else {
    loss <- lapply(lambdaSeq, valid_checkLoss, folds,
                   y=y, tau=tau, k=k)
  }
  
  loss <- matrix(unlist(loss), ncol =length(tau), byrow = T)
  total <- rowSums(loss)
  return(list(lambda = lambdaSeq[which.min(total)],
              loss = loss))
}
