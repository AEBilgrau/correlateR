context("rwishart and invwishart tests")

p <- 10
n <- 50
sigma <- diag(p)
sigma[rbind(1:2,2:1)] <- 0.5
nu <- 10


test_that("rwishart and rinvwishart returns proper arrays", {
  rw  <- rwishart(n, sigma, nu)
  expect_that(is.numeric(rw), is_true())
  expect_that(dim(rw), equals(c(p, p, n)))
  
  riw <- rinvwishart(n, sigma, nu)
  expect_that(is.numeric(riw), is_true())
  expect_that(dim(riw), equals(c(p, p, n)))
})




