context("Marginal covariance tests")

# Testing data
n <- 50
m_X <- 10
m_Y <- 15
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
  expect_equivalent(cov_XY,           xcovArma(X, Y, 0))
  expect_equivalent(cov_XY*(n - 1)/n, xcovArma(X, Y, 1))
  expect_equivalent(cov_XY,           xcovEigen(X, Y, 0))
  expect_equivalent(cov_XY*(n - 1)/n, xcovEigen(X, Y, 1))
  expect_equal(cov_XY,           xcov(X, Y))
  expect_equal(cov_XY*(n - 1)/n, xcov(X, Y, method = "ML"))
})

# Test for degenerate cases
combs <- expand.grid(n = c(0, 1, 10), m_X = c(0, 1, 10), m_Y = c(0, 1, 10))

for (i in seq_len(nrow(combs))) {
  n <- combs$n[i]
  m_X <- combs$m_X[i]
  m_Y <- combs$m_Y[i]
  X <- createData(n, m_X)
  Y <- createData(n, m_Y)
  cov_X <- stats::cov(X)
  cov_XY <- stats::cov(X, Y)

  test_that(sprintf("Covariance function works in degenerate cases, n = %d, m = %d", n, m_X), {
    expect_equivalent(cov_X, covArma(X))
    expect_equivalent(cov_X, covEigen(X))
    expect_equivalent(cov_X, covRcpp(X))
  })

  test_that(sprintf("Cross-covariance function works in degenerate cases, n = %d, m_X = %d, m_Y = %d", n, m_X, m_Y), {
    expect_equivalent(cov_XY, xcovArma(X, Y))
    expect_equivalent(cov_XY, xcovEigen(X, Y))
    expect_equivalent(cov_XY, xcovRcpp(X, Y))
  })
}
