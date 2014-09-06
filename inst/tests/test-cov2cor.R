context("cov2cor: conversion tests")

# Testing data
X    <- createData(n = 50, m = 10)
bigX <- createData(n = 5000, m = 1000)
S    <- cov(X)
bigS <- cov(bigX) 

ans <- stats::cov2cor(S)
res <- correlateR::cov2cor(S)

test_that("Correct conversion", {
  expect_that(ans,    equals(res))
  expect_that(cor(X), equals(res))
  expect_that(max(res) <=   1 + 2*.Machine$double.eps, is_true())
  expect_that(min(res) >= - 1 - 2*.Machine$double.eps, is_true())
})

test_that("Fast computation.", {
  expect_that(cov2cor(bigS),
              takes_less_than(system.time(stats::cov2cor(bigS))[3]))
})

# Add NA handling!!



