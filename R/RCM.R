#' Estimate degrees of freedom 
#' 
#' Function for estimating the degrees of freedom \eqn{nu}{\nu} in the 
#' random covariance model (RCM).
#' 
#' @param Psi A numeric matrix of size \eqn{p} times \eqn{p} giving the initial
#'   estimate of \eqn{Psi}{\Psi}.
#' @param S_list A \code{list} of scatter matrices.
#' @param ns Vector of group sizes.
#' @return A single number giving the \eqn{nu}{\nu} optimizing the RCM 
#'   likelihood with fixed \eqn{Psi}{\Psi}.
#' @author Anders Ellern Bilgrau
#' @keywords internal
rcm_get_nu <- function(Psi, S_list, ns) {
  loglik_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    rcm_loglik_nu_arma(Psi, nu, S_list, ns)
  }
  interval <- c(nrow(Psi) - 1 + 1e-10, 1e6) 
  res <- optimize(f = loglik_nu, interval = interval, maximum = TRUE)$maximum
  return(res)
} 

# Compute new Psi from nu, S, ns using pooled moment estimate
rcm_pool_step <- function(nu, S_list, ns, ...) {
  pool.sigma <- pool(S_list = S_list, ns = ns, norm_type = 1)
  fac <- nu - nrow(pool.sigma) - 1
  return(fac*pool.sigma)
}

# Compute new Psi from nu, S, ns using mean moment estimate
rcm_mean_step <- function(nu, S_list, ns, ...) {
  k <- length(ns)
  mean.sigma <-
    Reduce("+", lapply(seq_along(ns), function(i) S_list[[i]]/ns[i]))/k
  fac <- nu - nrow(mean.sigma) - 1
  return(fac*mean.sigma)
}

# Compute new Psi from nu, S, ns using approximate MLE
rcm_mle_step <- function(nu, S_list, ns, ...) {
  n.tot <- sum(ns)
  df <- nu + ns
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) df[i]*S_list[[i]]))/n.tot
  return(Psi)
}

#' Fit using the EM algorithm
#' 
#' Fit the RCM using the modified EM algorithm.
#' 
#' @param S A \code{list} of square scatter matrices of the same size.
#' @param ns A vector of group or sample sizes corresponding to the scatter 
#'   matrices in \code{S}.
#' @param Psi.init A \code{matrix} giving the initial estimate of 
#'   \eqn{Psi}{Psi}. Default starting value is the scaled pooled sample
#'   covariance matrix.
#' @param max.ite A numeric of length one giving the maximum number of 
#'   iterations allowed. Default is 1000.
#' @param nu.init A numeric of length one giving the inital estimate of 
#'   \eqn{nu}{nu}. Default is \code{sum(ns)}.
#' @param method The method to be used.
#' @param eps The convergence criterion.
#' @param verbose If true, the differences in log-likelihood for each iteration
#'   is printed out.
#' @return A named list of length 3 with the elements:
#'   \item{Psi}{A matrix giving the estimate of \eqn{Psi}{Psi}.}
#'   \item{nu}{A number giving the estimate of \eqn{nu}{nu}.}
#'   \item{iterations}{A integer giving the number of iterations used.}
#' @seealso \code{\link{Psi2Sigma}}
#' @examples
#' ns <- rep(11, 3)
#' print(Psi <- drop(rwishart(1)))
#' nu <- 30
#' S <- createRCMData(ns, Psi, nu)
#' 
#' with(fit.rcm(S, ns, method = "EM"),        Psi2Sigma(Psi, nu))
#' with(fit.rcm(S, ns, method = "pool"),      Psi2Sigma(Psi, nu))
#' with(fit.rcm(S, ns, method = "mean"),      Psi2Sigma(Psi, nu))
#' with(fit.rcm(S, ns, method = "approxMLE"), Psi2Sigma(Psi, nu))
#' Psi2Sigma(Psi, nu)
#' @export
fit.rcm <- function(S,
                    ns,
                    Psi.init,
                    nu.init,
                    method = c("EM", "pool", "mean", "approxMLE"),
                    max.ite = 1000, 
                    eps = 1e-3,
                    verbose = FALSE) {
  method <- match.arg(method)
  p <- nrow(S[[1]])
  if (missing(nu.init)) {
    nu.init <- sum(ns)
  }
  if (missing(Psi.init)) {
    Psi.init <- (nu.init - p - 1)*pool(S, ns)
  }
  updatePsi <- switch(method, 
                      "EM" = rcm_em_step_arma, "pool" = rcm_pool_step, 
                      "mean" = rcm_mean_step, "approxMLE" = rcm_mle_step)
  if (method == "EM") {
    conv <- function(x) {
      if (x < 0) warning("log-likelihood increased in last iteration!")
      return(x)
    }
  } else {
    conv <- abs
  }
  
  Psi.old <- Psi.init
  nu.old <- nu.init
  ll.old  <- rcm_loglik_arma(Psi = Psi.old, nu = nu.old, S_list = S, ns = ns)
  for (i in seq_len(max.ite)) {
    Psi.new <- updatePsi(Psi = Psi.old, nu = nu.old, S_list = S, ns = ns)
    nu.new  <- rcm_get_nu(Psi = Psi.new, S_list = S, ns = ns)
    ll.new  <- rcm_loglik_arma(Psi = Psi.new, nu = nu.new, S_list = S, ns = ns)
    if (conv(ll.new - ll.old) < eps) {
      break
    } else {
      if (verbose) {
        cat(sprintf("it = %g : L = %.3f : diff = %.3f\n", 
                    i, ll.new, ll.new - ll.old))
      }
      Psi.old <- Psi.new
      nu.old  <- nu.new
      ll.old  <- ll.new
    }
  }
  if (i == max.ite) warning("max iterations (", max.ite, ") hit!")
  return(list("Psi" = Psi.new, "nu" = nu.new, "iterations" = i))
}

#' Conversion from Psi and nu to Sigma
#' 
#' Computes the expected covariance matrix from Psi and nu in the random 
#' covariance model (RCM).
#' 
#' @param Psi A numeric square positive semi-definite matrix.
#' @param nu A numeric of length one giving the degrees of freedom in the RCM.
#' @return A covariance matrix the same size as \code{Psi}.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @export
Psi2Sigma <- function(Psi, nu) {
  return(Psi/(nu - ncol(Psi) - 1))
}

#' Density of the RCM model
#' 
#' The density of the RCM model when the mean \code{mu} is non-zero. Used for
#' the HDA analysis.
#' 
#' @param x A numeric matrix of size n by m with observations in the rows and
#'   variables in the columns.
#' @param mu A numeric vector of length m giving the mean values of the RCM in
#'   HDA.
#' @param Psi A numeric matrix giving the parameter of the RCM.
#' @param nu A numeric of length 1 giving the degrees of freedom.
#' @param logarithm A boolean of length 1. If \code{TRUE}, the log density is 
#'   returned. Defaults to \code{FALSE}.
#' @return 
#'   Returns a vector of length n where the \eqn{i}th corresponds to the density
#'   evaluated in the \eqn{i}th row in \code{x}.
#' @examples
#' n <- 10
#' m <- 5
#' x <- matrix(rnorm(n*m), n, m)
#' mu <- rnorm(m)
#' Psi <- cov(matrix(rnorm(n*m), n, m))
#' nu <- 15
#' drcm(x, mu, Psi, nu)
#' drcm(x, mu, Psi, nu, logarithm = TRUE)
#' @export
drcm <- function(x, mu, Psi, nu, logarithm = FALSE) {
  p <- nrow(Psi)
  if (is.null(dim(x))) {
    dim(x) <- c(1, p)
  }
  Q <- function(x, A) {
    rowSums(tcrossprod(x, A) * x)
  }
  stopifnot(ncol(x) == length(mu))
  t1 <- lgammap((nu + 1)/2, p = p)
  t2 <- -lgammap(nu/2, p = p)
  t3 <- p/2*log(pi)
  t4 <- -1/2*logdet_arma(Psi)[1]
  
  x.center <- t(t(x) - mu)
  t5 <- -(nu + 1)/2 * log(1 + Q(x.center, solve(Psi)))
  
  ans <- t1 + t2 + t3 + t4 + t5
  if (!logarithm) {
    ans <- exp(ans)
  }
  attributes(ans) <- NULL
  return(ans)
}

#' Simulate data from RCM
#' 
#' Generate data from the hierarchical random covariance model (RCM).
#' 
#' @param ns A numeric vector giving the sample sizes in each study.
#' @param psi The underlying \eqn{Psi} parameter. If \code{nu} is \code{Inf} 
#'   this is used as \eqn{Sigma} parameter in all Wishart distribution.
#' @param nu A numeric of length one giving the underlying \eqn{nu} parameter.
#' @return A \code{list} of matrices of the same size as \code{psi} giving
#'   observed scatter matrices from the RCM.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @examples
#' ns <- c(20, 14, 10)
#' psi <- diag(3)
#' createRCMData(ns, psi, nu = 10)
#' 
#' createRCMData(ns, psi, nu = 1e20)
#' # is NOT the same as
#' createRCMData(ns, psi, nu = Inf)
#' # which is almost the same as
#' createRCMData(ns, psi*(1e20 - 3 - 1), nu = 1e20)
#' @export
createRCMData <- function(ns, psi, nu) {
  stopifnot(length(ns) > 0)
  k <- length(ns)
  if (nu == Inf) {
    S <- lapply(seq_len(k), function(i) drop(rwishart(1, psi, ns[i])))
  } else {
    sigmas <- rinvwishart(n = k, psi = psi, nu = nu)
    S <- lapply(seq_len(k), function(i) drop(rwishart(1, sigmas[, , i], ns[i])))
  }
  for (i in seq_along(S)) {
    dimnames(S[[i]]) <- dimnames(psi)
  }
  return(S)
}


