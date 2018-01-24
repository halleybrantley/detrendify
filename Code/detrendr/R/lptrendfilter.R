lptrendfilter <- function(y, tau, lambda, k){
  n <- length(y)
  D <- get_Dk(length(y), k)
  m <- nrow(D)
  f.obj <- c(rep(lambda, 2*m), rep(tau, n), rep(1-tau, n))  
  f.con <- as.matrix(rbind(cbind(diag(m),matrix(0, nrow=m, ncol=m+2*n)),
                           cbind(matrix(0, nrow=m, ncol=m), diag(m), 
                                 matrix(0, nrow=m, ncol=2*n)),
                           cbind(matrix(0, nrow=n, ncol=2*m), diag(n), 
                                 matrix(0, nrow=n, ncol=n)), 
                           cbind(matrix(0, nrow=n, ncol=2*m+n), diag(n)),
                           cbind(diag(m), -diag(m), D, -D)))
  f.dir <- c(rep(">=", nrow(f.con)-m), rep("=", m))
  f.rhs <- c(rep(0, nrow(f.con)-m), as.numeric(D%*%y))
  lp.fit <- lp("min", f.obj, f.con, f.dir, f.rhs)
  gamma <- lp.fit$solution
  theta <- y - gamma[(2*m+1):(2*m+n)] + gamma[(2*m+n+1):(2*(m+n))]
  return(theta)
}