
test <- function(X, ...) {
  UseMethod("test")
}

a <- 1
attr(a, "class") <- "cor"

test.cor <- function(X) {
  cat("test.cor called\n") 
  return(-1)
}

test(a)

