context("Marginal correlation tests")

# Testing data
n <- 50
X <- replicate(10, rnorm(n))
dimnames(X) <- list(paste0("obs", 1:n), paste0("dim", 1:ncol(X)))
Y <- replicate(15, rnorm(n))
dimnames(Y) <- list(paste0("obs", 1:n), paste0("dim", 1:ncol(Y)))

# Gold standard
cor_X  <- stats::cor(X)     # Standard (marginal) correlation
cor_XY <- stats::cor(X, Y)  # Standard (marginal) cross-correlation

test_that("Correct correlation.", {
  expect_that(cor_X, is_equivalent_to( corArma(X)))
  expect_that(cor_X, is_equivalent_to(xcorArma(X, X)))
  expect_that(cor_X, is_equivalent_to( corEigen(X)))
  expect_that(cor_X, is_equivalent_to(xcorEigen(X, X)))
  expect_that(cor_X, is_equivalent_to( corRcpp(X)))
  expect_that(cor_X, is_equivalent_to(xcorRcpp(X, X)))
  expect_that(cor_X, equals(cor(X)))
})

test_that("Correct cross-correlation.", {
  expect_that(cor_XY, is_equivalent_to(xcorArma(X, Y)))
  expect_that(cor_XY, is_equivalent_to(xcorEigen(X, Y)))
  expect_that(cor_XY, is_equivalent_to(xcorRcpp(X, Y)))
  expect_that(cor_XY, equals(xcor(X, Y)))
})


#
# With NAs
#


# Add random NAs
X[sample(length(X), 10)] <- NA
Y[sample(length(Y), 10)] <- NA

# Recompute correlation standard
# Note, unlike cor(X), cor(X, X) returns NAs in the diagonal if NAs are present
# in X. correlatoR also does this! 
cor_X  <- stats::cor(X, X) # Standard (marginal) correlation 
cor_XY <- stats::cor(X, Y) # Standard (marginal) cross-correlation

test_that("Correct correlation with NAs.", {
  expect_that(cor_X, is_equivalent_to( corArma(X)))
  expect_that(cor_X, is_equivalent_to(xcorArma(X, X)))
  expect_that(cor_X, is_equivalent_to( corEigen(X)))
  expect_that(cor_X, is_equivalent_to(xcorEigen(X, X)))
  expect_that(cor_X, is_equivalent_to( corRcpp(X)))
  expect_that(cor_X, is_equivalent_to(xcorRcpp(X, X)))
  expect_that(cor_X, equals(cor(X)))
})

test_that("Correct cross-correlation with NAs.", {
  expect_that(cor_XY, is_equivalent_to(xcorArma(X, Y)))
  expect_that(cor_XY, is_equivalent_to(xcorEigen(X, Y)))
  expect_that(cor_XY, is_equivalent_to(xcorRcpp(X, Y)))
  expect_that(cor_XY, equals(xcor(X, Y)))
})



#
# Tests for degenerate cases
#

combs <- expand.grid(n = c(0, 1, 10), m_X = c(0, 1, 10), m_Y = c(0, 1, 10))

for (i in seq_len(nrow(combs))) {
  n <- combs$n[i]
  m_X <- combs$m_X[i]
  m_Y <- combs$m_Y[i]
  X <- createData(n, m_X)
  Y <- createData(n, m_Y)
  cor_X <- stats::cor(X)
  cor_XY <- stats::cor(X, Y)
  
  test_that(sprintf("corariance function works in degenerate cases, n = %d, m = %d", n, m_X), {
    expect_equivalent(cor_X, corArma(X))
    expect_equivalent(cor_X, corEigen(X))
    expect_equivalent(cor_X, corRcpp(X))
  })
  
  test_that(sprintf("Cross-corariance function works in degenerate cases, n = %d, m_X = %d, m_Y = %d", n, m_X, m_Y), {
    expect_equivalent(cor_XY, xcorArma(X, Y))
    expect_equivalent(cor_XY, xcorEigen(X, Y))
    expect_equivalent(cor_XY, xcorRcpp(X, Y))
  })
}
