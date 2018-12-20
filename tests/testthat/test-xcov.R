context("Marginal covariance tests")

# Testing data
n <- 52
m_X <- 11
m_Y <- 13
X <- createData(n, m_X)
Y <- createData(n, m_Y)

# Gold standard:
xcov_X  <- stats::cov(X)  # Standard covariance
xcov_XY <- stats::cov(X, Y)  # Standard cross-covariance

test_that("Dimensions are as expected", {
  expect_equal(dim(xcov_X),  c(m_X, m_X))
  expect_equal(dim(xcov_XY), c(m_X, m_Y))
})


test_that("Covariance is computed correctly", {
  expect_equivalent(xcov_X,           covArma(X, 0))
  expect_equivalent(xcov_X*(n - 1)/n, covArma(X, 1))
  expect_equivalent(xcov_X,           xcovArma(X, X, 0))
  expect_equivalent(xcov_X*(n - 1)/n, xcovArma(X, X, 1))
  expect_equivalent(xcov_X,           covEigen(X, 0))
  expect_equivalent(xcov_X*(n - 1)/n, covEigen(X, 1))
  expect_equivalent(xcov_X,           xcovEigen(X, X, 0))
  expect_equivalent(xcov_X*(n - 1)/n, xcovEigen(X, X, 1))
  expect_equal(xcov_X,           cov(X))
  expect_equal(xcov_X*(n - 1)/n, cov(X, method = "ML"))
})

test_that("Cross-correlation is computed correctly", {
  expect_equivalent(xcov_XY,           xcovArma(X, Y, 0))
  expect_equivalent(xcov_XY*(n - 1)/n, xcovArma(X, Y, 1))
  expect_equivalent(xcov_XY,           xcovEigen(X, Y, 0))
  expect_equivalent(xcov_XY*(n - 1)/n, xcovEigen(X, Y, 1))
  expect_equivalent(xcov_XY,           xcov(X, Y))
  expect_equivalent(xcov_XY*(n - 1)/n, xcov(X, Y, method = "ML"))
})

# ADD TESTS WITH MISSING VALUES!



