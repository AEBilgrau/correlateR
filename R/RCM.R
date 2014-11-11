#' Estimate degrees of freedom 
#' 
#' Function for estimating the degrees of freedom \eqn{nu}{\nu} in the 
#' random covariance model (RCM).
#' 
#' @param Psi A numeric matrix of size \eqn{p} times \eqn{p} giving the initial
#'   estimate of \eqn{Psi}{\Psi}.
#' @param nu  A single numeric number giving the degrees of freedom 
#'   \eqn{nu}{\nu}.
#' @param S A \code{list} of scatter matrices.
#' @param ns Vector of group sizes.
#' @return A single number giving the \eqn{nu}{\nu} optimizing the RCM 
#'   likelihood with fixed \eqn{Psi}{\Psi}.
#' @author Anders Ellern Bilgrau
#' @keywords internal
rcm_get_nu <- function(Psi, nu, S, ns, interval) {
  # Find maxima with optimize
  loglik_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    rcm_loglik_nu_arma(Psi, nu, S, ns)
  }
  interval <- c(nrow(Psi) - 1 + 1e-10, 1e6)
  res <- optimize(f = loglik_nu, interval = interval, maximum = TRUE)$maximum
  return(res)
} 

# Compute new Psi from nu, S, ns using moment estimate
rcm_moment_step <- function(nu, S, ns) {
  k <- length(ns)
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) S[[i]]/ns[i]))/k
  p <- nrow(Psi)
  fac <- nu - p - 1
  return(fac*Psi)
}

# Compute new Psi from nu, S, ns using approximate MLE
rcm_mle_step <- function(nu, S, ns) {
  n.tot <- sum(ns)
  fac <- nu + ns
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) fac[i]*S[[i]]))/n.tot
  return(Psi)
}

#' Fit using the EM algorithm
#' 
#' Fit the RCM using the modified EM algorithm.
#' 
#' @param S A \code{list} of square scatter matrices of the same size.
#' @param ns A vector of sample sizes corresponding to the scatter matrices in 
#'   \code{S}.
#' @param Psi.init A \code{matrix} giving the initial estimate of 
#'   \eqn{Psi}{Psi}.
#' @param max.ite A numeric of length one giving the maximum number of 
#'   iterations allowed.
#' @param nu.init A numeric of length one giving the inital estiamte of 
#'   \eqn{nu}{nu}.
#' @param eps The convergence criterion.
#' @param verbose If true, the differences in log-likelihood for each iteration
#'   is printed out.
#' @return A named list of length 3 with the elements:
#'   \item{Psi}{A matrix giving the estimate of \eqn{Psi}{Psi}.}
#'   \item{nu}{A number giving the estimate of \eqn{nu}{nu}.}
#'   \item{iterations}{A integer giving the number of iterations used.}
#' @seealso \code{\link{Psi2Sigma}}
#' @export
fit.rcm <- function(S,
                     ns,
                     Psi.init = correlateR:::pool(S, ns),
                     nu.init = sum(ns) + 1,
                     max.ite = 1000, 
                     eps = 1e-3,
                     verbose = FALSE) {
  p <- nrow(S)
  interval <- c(p - 1 + 10*.Machine$double.eps, 1e6)
  Psi.old <- Psi.init
  nu.old  <- nu.init
  for (i in seq_len(max.ite)) {
    ll.old  <- rcm_loglik_arma(Psi.old, nu.old, S, ns)
    Psi.new <- rcm_em_step_arma(Psi.old, nu.old, S, ns)
    nu.new  <- rcm_get_nu(Psi.new, nu.old, S, ns, interval)
    ll.new  <- rcm_loglik_arma(Psi.new, nu.new, S, ns)
    stopifnot(ll.new > ll.old)
    if (ll.new - ll.old < eps) {
      break
    } else {
      Psi.old <- Psi.new
      nu.old <- nu.new
    }
    if (verbose) {
      cat("ite =", i, ":", "ll.new - ll.old =", ll.new - ll.old, "\n");
      flush.console()
    }
  }
  if (i == max.ite) warning("max iterations (", max.ite, ") hit!")
  return(list("Psi" = Psi.new, "nu" = nu.new, "iterations" = i))
}


# MLE alg
fit.rcm.MLE <- function(S, ns,
                         nu.init = nrow(S[[1]]) + 2,
                         max.ite = 1000, eps = 1e-3,
                         verbose = FALSE) {
  p <- nrow(S)
  interval <- c(p - 1 + 10*.Machine$double.eps, 1e6)
  nu.old <- nu.init
  Psi.old <- rcm_mle_step(nu = nu.old, S = S, ns = ns)
  for (i in seq_len(max.ite)) {
    ll.old <- rcm_loglik_arma(Psi.old, nu.old, S, ns)
    nu.new  <- rcm_get_nu(Psi.old, nu.old, S, ns, interval)
    Psi.new <- rcm_mle_step(nu.new, S, ns)
    ll.new  <- rcm_loglik_arma(Psi.new, nu.new, S, ns)
    stopifnot(ll.new > ll.old)
    if (ll.new - ll.old < eps) {
      break
    } else {
      Psi.old <- Psi.new
      nu.old  <- nu.new
    }
    if (verbose) {
      cat("ite =", i, ":", "ll.new - ll.old =", ll.new - ll.old, "\n");
      flush.console()
    } 
  }
  
  if (i == max.ite) warning("max iterations (", max.ite, ") hit!")
  
  return(list("Psi" = Psi.new, "nu" = nu.new, "iterations" = i))
}

# Estimation using moment
fit.rcm.moment <- function(S, ns,
                            nu.init = nrow(S[[1]]) + 10,
                            max.ite = 1000, eps = 1e-3,
                            verbose = FALSE) {
  p <- nrow(S[[1]])
  interval <- c(p - 1 + 10*.Machine$double.eps, 1e6)
  nu.old   <- nu.init
  Psi.old  <- rcm_moment_step(nu = nu.old, S = S, ns = ns)
  for (i in seq_len(max.ite)) {
    nu.new  <- rcm_get_nu(Psi = Psi.old, nu = nu.old, S = S, ns = ns,
                           interval = interval)
#     fac <- (nu.new - p - 1)/(nu.old - p - 1)
#     Psi.new <- Psi.old * fac
    Psi.new  <- rcm_moment_step(nu = nu.old, S = S, ns = ns)
    if (verbose) {
      cat("ite =", i, ":", "nu.new - nu.old =", nu.new - nu.old,
          "nu =", nu.new, "\n");
      flush.console()
    }
    if (abs(nu.new - nu.old) < eps) {
      break
    } else {
      Psi.old <- Psi.new
      nu.old  <- nu.new
    }
  }
  if (i == max.ite) warning("max iterations (", max.ite, ") hit!")
  return(list("Psi" = Psi.new, "nu" = nu.new, "iterations" = i))
}

# Conversion from Psi and nu to Sigma
 
#' @export
Psi2Sigma <- function(Psi, nu) {
  return(Psi/(nu - ncol(Psi) - 1))
}

# Density of the RCM model
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
