context("Projection onto V")

test_that("Projection onto V produces expected result", {
  set.seed(12345)
  k <- 3
  n <- 1e2
  D <- get_Dk_R(n,k)
  theta <- rnorm(n)
  eta <- as.numeric(D%*%theta)
  L1 <- test_project_V(theta, eta, n, k)
  
  test_theta(test_project_V(theta, eta, n, k), 
             equivalent_to())
  
})