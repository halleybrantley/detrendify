phiBar_list <- list()
for (i in 1:n_windows){
  phiBar_list[[i]] <- as.numeric(phiBar[windows[,i],])
}

# Window Update
update_loop <- function(y_list, w_list, phiBar_list, tau, lambda, k, rho){
  phi_list <- vector("list",length(y_list))
  
  for (i in 1:n_windows){
    phi_list[[i]] <- quad_update(y_list[[i]], tau, lambda, k, w_list[[i]],
                                 phiBar_list[[i]], rho)
  }
  return(phi_list)
}

update_mapply <- function(y_list, w_list, phiBar_list, tau, lambda, k, rho){
  
  phi_list <- mapply(quad_update, y_list, w_list, phiBar_list,
                     MoreArgs = list(tau=tau, lambda=lambda, k=k,
                                     rho = rho),
                     SIMPLIFY = FALSE)
  return(phi_list)
}

library(parallel)
update_mcapply <- function(y_list, w_list, phiBar_list, tau, lambda, k, rho){
  
  phi_list <- mcmapply(quad_update, y_list, w_list, phiBar_list,
                       MoreArgs = list(tau=tau, lambda=lambda, k=k,
                                       rho = rho),
                       SIMPLIFY = FALSE)
  return(phi_list)
}



microbenchmark(
  #update_loop(y_list, w_list, phiBar_list, tau, lambda, k, rho),
  update_mapply(y_list, w_list, phiBar_list, tau, lambda, k, rho),
  update_mcapply(y_list, w_list, phiBar_list, tau, lambda, k, rho),
  update_foreach(y_list, w_list, phiBar_list, tau, lambda, k, rho),
  times = 5
)

stopCluster(cl)
