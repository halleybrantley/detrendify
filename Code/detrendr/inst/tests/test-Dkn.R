context("Discrete Derivatives")

test_that("Dk(n) produces expected output",{
  n <- 10
  k <- 4
  expect_that(test_Dk(n, k, 0, 0), equals(get_Dk_R(n,k)[1,1]))
  expect_that(test_Dk(n, k, 1, 3), equals(get_Dk_R(n,k)[2,4]))
  expect_that(test_Dk(n, k, 2, 4), equals(get_Dk_R(n,k)[3,5]))
})
