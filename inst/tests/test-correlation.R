context("Marginal correlation tests")

# Testing data
n <- 50
X <- replicate(10, rnorm(n))
dimnames(X) <- list(paste0("obs", 1:n), paste0("dim", 1:ncol(X)))
Y <- replicate(15, rnorm(n))
dimnames(Y) <- list(paste0("obs", 1:n), paste0("dim", 1:ncol(Y)))

# Gold standard:
cor_X <- stats::cor(X)  # Standard (marginal) correlation
cor_XY <- stats::cor(X, Y)  # Standard (marginal cross-correlation
all.equal(cor_X*(n - 1)/n, corEigen(X, 1))


test_that("Correct computation of correlation.", {
  expect_that(cor_X,           is_equivalent_to(corArma(X, 0)))
  expect_that(cor_X*(n - 1)/n, is_equivalent_to(corArma(X, 1)))
  expect_that(cor_X,           is_equivalent_to(crosscorArma(X, X, 0)))
  expect_that(cor_X*(n - 1)/n, is_equivalent_to(crosscorArma(X, X, 1)))
  expect_that(cor_X,           is_equivalent_to(corEigen(X, 0)))
  expect_that(cor_X*(n - 1)/n, is_equivalent_to(corEigen(X, 1)))
  expect_that(cor_X,           is_equivalent_to(crosscorEigen(X, X, 0)))
  expect_that(cor_X*(n - 1)/n, is_equivalent_to(crosscorEigen(X, X, 1)))
  expect_that(cor_X,           equals(corariance(X)))
  expect_that(cor_X*(n - 1)/n, equals(corariance(X, method = "ML")))
})

test_that("Correct computation of cross-correlation.", {
  expect_that(cor_XY,           is_equivalent_to(crosscorArma(X, Y, 0)))
  expect_that(cor_XY*(n - 1)/n, is_equivalent_to(crosscorArma(X, Y, 1)))
  expect_that(cor_XY,           is_equivalent_to(crosscorEigen(X, Y, 0)))
  expect_that(cor_XY*(n - 1)/n, is_equivalent_to(crosscorEigen(X, Y, 1)))
  expect_that(cor_XY,           equals(corariance(X, Y)))
  expect_that(cor_XY*(n - 1)/n, equals(corariance(X, Y, method = "ML")))
})

