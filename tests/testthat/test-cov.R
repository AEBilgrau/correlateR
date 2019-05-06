context("Marginal covariance tests")

# Testing data
n <- 50
m_X <- 10
m_Y <-15
X <- createData(n, m_X)
Y <- createData(n, m_Y)

# Gold standard:
cov_X <- stats::cov(X)  # Standard covariance
cov_XY <- stats::cov(X, Y)  # Standard cross-covariance

test_that("Covariance is computed correctly", {
  expect_that(cov_X,           is_equivalent_to(covArma(X, 0)))
  expect_that(cov_X*(n - 1)/n, is_equivalent_to(covArma(X, 1)))
  expect_that(cov_X,           is_equivalent_to(xcovArma(X, X, 0)))
  expect_that(cov_X*(n - 1)/n, is_equivalent_to(xcovArma(X, X, 1)))
  expect_that(cov_X,           is_equivalent_to(covEigen(X, 0)))
  expect_that(cov_X*(n - 1)/n, is_equivalent_to(covEigen(X, 1)))
  expect_that(cov_X,           is_equivalent_to(xcovEigen(X, X, 0)))
  expect_that(cov_X*(n - 1)/n, is_equivalent_to(xcovEigen(X, X, 1)))
  expect_that(cov_X,           equals(cov(X)))
  expect_that(cov_X*(n - 1)/n, equals(cov(X, method = "ML")))
})

test_that("Cross-correlation is computed correctly", {
  expect_that(cov_XY,           is_equivalent_to(xcovArma(X, Y, 0)))
  expect_that(cov_XY*(n - 1)/n, is_equivalent_to(xcovArma(X, Y, 1)))
  expect_that(cov_XY,           is_equivalent_to(xcovEigen(X, Y, 0)))
  expect_that(cov_XY*(n - 1)/n, is_equivalent_to(xcovEigen(X, Y, 1)))
  expect_that(cov_XY,           equals(xcov(X, Y)))
  expect_that(cov_XY*(n - 1)/n, equals(xcov(X, Y, method = "ML")))
})

# Test for degenerate cases
combs <- expand.grid(n = c(0, 1, 10), m_X = c(0, 1, 10), m_Y = c(0, 1, 10))

for (i in seq_len(nrow(combs))) {
  n <- combs$n[i]
  m_X <- combs$m_X[i]
  m_Y <- combs$m_Y[i]
  x <- createData(n, m_X)
  y <- createData(n, m_Y)
  cov_x <- stats::cov(x)
  cov_xy <- stats::cov(x, y)
  
  test_that(sprintf("Covariance function works in degenerate cases, n = %d, m = %d", n, m_X), {
    expect_equal(cov_x, covArma(x))
    expect_equal(cov_x, covArma(x))
    expect_equal(cov_x, covRcpp(x)) 
  })
  
  test_that(sprintf("Cross-covariance function works in degenerate cases, n = %d, m_X = %d, m_Y = %d", n, m_X, m_Y), {
    expect_equal(cov_xy, covArma(x))
    expect_equal(cov_xy, covArma(x))
    expect_equal(cov_xy, covRcpp(x)) 
  })
}
