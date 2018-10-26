#' Estimate trend using rolling windows
#'
#' \code{get_windows} Estimates trend in rolling windows
#'
#' @param y observed data
#' @param time timestamp for observations
#' @param k order of differencing
#' @param tau quantile levels at which to evaluate trend
#' @param lambda smoothing penalty parameter
#' @param window_size integer length of window
#' @param overlap integer length of overlap between windows
get_windows <- function(y, time, k, tau, lambda, 
                        window_size, overlap){
  y_n <- length(y)
  window_end <- window_size
  window_start <- 1
  
  theta_df <- {}
  i <- 1
  pb <- txtProgressBar(min = 0, max = round(y_n/(window_size-overlap)),
                       style = 3)
  
  while(window_end < y_n){
    y0 <- y[window_start:window_end]
    theta_gurobi <- as.data.frame(gurobi_trend(y0, tau, lambda, k))
    names(theta_gurobi) <- paste0("tau_", tau)
    theta_gurobi$time <- time[window_start:window_end]
    theta <- theta_gurobi %>% gather(tau, theta, -time)
    theta$window <- i
    theta_df <- bind_rows(theta_df, theta)
    window_start <- window_start + window_size - overlap
    window_end <- window_start + window_size - 1
    i <- i+1
    setTxtProgressBar(pb, i)
  }
  
  window_end <- y_n
  y0 <- y[window_start:window_end]
  theta_gurobi <- as.data.frame(gurobi_trend(y0, tau, lambda, k))
  names(theta_gurobi) <- paste0("tau_", tau)
  theta_gurobi$time <- time[window_start:window_end]
  theta <- theta_gurobi %>% gather(tau, theta, -time)
  theta$window <- i
  theta_df <- bind_rows(theta_df, theta)
  
  return(theta_df)
}
