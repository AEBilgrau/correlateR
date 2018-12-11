context("rwishart and invwishart tests")

p <- 10
n <- 50
sigma <- diag(p)
sigma[rbind(1:2,2:1)] <- 0.5
nu <- 10

rw  <- rwishart(n, sigma, nu)
riw <- rinvwishart(n, sigma, nu)

test_that("rwishart and rinvwishart returns proper array", {
  expect_that(is.numeric(rw), is_true())
  expect_that(dim(rw), equals(c(p, p, n)))
  expect_that(is.numeric(riw), is_true())
  expect_that(dim(riw), equals(c(p, p, n)))
})




