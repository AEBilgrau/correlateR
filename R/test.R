
test <- function(X, ...) {
  UseMethod("test")
}

test.cor <- function(X) {
  cat("Test.cor called\n") 
  return(1)
}


# a <- matrix(1, 2, 2)
# str(a)
# class(a)
# 
# attr(a, "class") <- "cor"
# inherits(a, "matrix")
# 
# class(a)
# Test(a)

