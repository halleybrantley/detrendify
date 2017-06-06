context("Projection onto V")

test_that("Projection onto V produces expected result", {
  set.seed(12345)
  k <- 3
  n <- 1e2
  D <- get_Dk_R(n,k)
  M <- diag(n) + crossprod(D)
  theta <- rnorm(n)
  eta <- as.numeric(D%*%theta) + 0.01*rnorm(n-k)
  #L1 <- test_project_V(theta, eta, n, k)
  #L2 <- project_V_R(theta, eta, D, M)
  expect_that(test_project_V(theta, eta, k), 
             is_equivalent_to(project_V_R(theta, eta, D, M)))
  
})

