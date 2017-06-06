context("Spingarn Step")

test_that("Single Spingarn step produces expected result", {
  set.seed(12345)
  n <- 1e2
  x <- seq(1/n, 1, length.out=n)
  f <- 2*(x + 2)^2 + 3*cos(3*pi*x)
  tau <- 1e4
  g <-100*exp(-tau*(x-0.5)^2)
  y <- f + g + rnorm(n)
  k <- 3
  D <- get_Dk_R(n, k)
  M <- diag(n) + crossprod(D)
  lambda <- 1
  tau <- 0.01
  step <- 1
  theta <- y
  eta <- as.numeric(matrix(D %*% theta))
  
  expect_that(spingarn_one_step_R(theta, eta, y, D, lambda, tau, step, M), 
              is_equivalent_to(
                spingarn_multi_iter(theta, eta, y, n, k, lambda, tau, step, 1)))
})
  
