#' @rdname corcov
#' @export
pxcor <- function(X, Y, Z) {
  pxcorSimple <- function(x, y, Z) {
    Z1 <- cbind(1, Z)
    r_xz <- fastLmPure(Z1, x)$residuals # Alternative: residuals(lm(x ~ Z))
    r_yz <- fastLmPure(Z1, y)$residuals
    return(cor(r_xz, r_yz))
  }
  
  ans <- matrix(NA, ncol(X), ncol(Y))
  for (i in 1:ncol(X)) {
    for (j in 1:ncol(Y)) {
      ans[i,j] <- pxcorSimple(X[, i], Y[, j], Z = Z)
    }
  }
  return(ans)
}


