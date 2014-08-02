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
  expect_that(cor_X, is_equivalent_to(corArma(X)))
  expect_that(cor_X, is_equivalent_to(crosscorArma(X, X)))
  expect_that(cor_X, is_equivalent_to(corEigen(X)))
  expect_that(cor_X, is_equivalent_to(crosscorEigen(X, X)))
  expect_that(cor_X, is_equivalent_to(corRcpp(X)))
  expect_that(cor_X, is_equivalent_to(crosscorRcpp(X, X)))
  expect_that(cor_X, equals(correlation(X)))
})

test_that("Correct cross-correlation.", {
  expect_that(cor_XY, is_equivalent_to(crosscorArma(X, Y)))
  expect_that(cor_XY, is_equivalent_to(crosscorEigen(X, Y)))
  expect_that(cor_XY, is_equivalent_to(crosscorRcpp(X, Y)))
  expect_that(cor_XY, equals(correlation(X, Y)))
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
cor_X  <- stats::cor(X, X)   # Standard (marginal) correlation 
cor_XY <- stats::cor(X, Y)  # Standard (marginal) cross-correlation

test_that("Correct correlation with NAs.", {
  expect_that(cor_X, is_equivalent_to(corArma(X)))
  expect_that(cor_X, is_equivalent_to(crosscorArma(X, X)))
  expect_that(cor_X, is_equivalent_to(corEigen(X)))
  expect_that(cor_X, is_equivalent_to(crosscorEigen(X, X)))
  expect_that(cor_X, is_equivalent_to(corRcpp(X)))
  expect_that(cor_X, is_equivalent_to(crosscorRcpp(X, X)))
  expect_that(cor_X, equals(correlation(X)))
})

test_that("Correct cross- correlation with NAs.", {
  expect_that(cor_XY,           is_equivalent_to(crosscorArma(X, Y)))
  expect_that(cor_XY,           is_equivalent_to(crosscorEigen(X, Y)))
  expect_that(cor_XY,           is_equivalent_to(crosscorRcpp(X, Y)))
  expect_that(cor_XY,           equals(correlation(X, Y)))
})
