#' Compute first order partial correlations
#'
#' This function computes the first order partial correlation given the
#' \code{k}'th variable.
#' 
#' @param S A numeric correlation matrix
#' @param k An integer. The index of the variable to condition on.
#' @return A matrix of the same size as \code{X} with \code{NaN}s in the
#'   \code{k} row and column.
#' @seealso \code{\link{cor}}
#' @author Anders Ellernsta Bilgrau <abilgrau (at) math.aau.dk>
#' @examples
#' X <- replicate(4, rnorm(10))
#' dimnames(X) <- list(paste0("obs", 1:nrow(X)), paste0("var", 1:ncol(X)))
#' S <- cor(X)
#' pcor(S, 2)
firstOrderPartialCorrelation <- function(S, k) {
  stopifnot(k %in% seq_len(ncol(S)))
  
  # 1 - S[k, ]^2 sometimes is negative!
  tmp <- sqrt(abs(1 - S[k, ]^2))
  
  ans <- (S - tcrossprod(S[k, ]))/tcrossprod(tmp)
  ans[k, ] <- ans[, k] <- NaN
  return(ans)
}


firstOrderPartialCorrelation2 <- function(X, Y, Z) {

  xcor <- cor(X, Y)
  xcor
  
  # 1 - S[k, ]^2 sometimes is negative!
  tmp <- sqrt(abs(1 - S[k, ]^2))
  
  ans <- (S - tcrossprod(S[k, ]))/tcrossprod(tmp)
  ans[k, ] <- ans[, k] <- NaN
  return(ans)
}
