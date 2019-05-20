context("Test centering matrix cen(d)")

for (d in 0:10) {
  X <- matrix(runif(d*d), d, d)
  test_that(sprintf("Checking cen(%d) works as intended", d), {
    cc <- cen(d)
    expect_true(is.matrix(cc))
    expect_equal(dim(cc), c(d, d))
    expect_equal(cc, t(cc)) # Symmetric
    expect_equal(crossprod(cc), cc) # Idempotent
    
    # Check for expected error
    expect_error(cen(d + runif(1))) # Non integer input
    expect_error(cen(d:30)) # Vector input
    
    # That the centering matrix does its job
    expect_equal(colMeans(cc %*% X), rep(0, d))
    expect_equal(rowMeans(X %*% cc), rep(0, d))
  })
  
}