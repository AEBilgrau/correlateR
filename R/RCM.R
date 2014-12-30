#' @rdname RCMmisc
#' @title RCM miscellaneous functions
#' @description Miscellaneous functions for the random covariance model (RCM). 
#' @details \code{ICC} compute the ICC in the RCM.
#' A simple function for computing the intra-class correlation coefficient 
#' (ICC). This function simply computes 1 divided by \code{nu - p}.
#' @param nu A numeric of length one giving the degrees of freedom in the RCM.
#' @param p A numeric giving the dimension of the space.
#' @return \code{ICC}: A numeric giving the ICC.
#' @export
ICC <- function(nu, p) {
  return(1/(nu - p))
}

#' @rdname RCMmisc
#' @details \code{Psi2Sigma} and \code{Sigma2Psi} provide conversion between Psi 
#' and Sigma.
#' Computes the expected covariance matrix from Psi and nu in the random 
#' covariance model (RCM) and vice versa.
#' @param Psi A numeric square positive semi-definite matrix. The underlying 
#'   parameter in the RCM.
#' @return \code{Psi2Sigma}, \code{Sigma2Psi}:
#'   The converted matrix the same size as \code{Psi} or \code{Sigma}.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @export
Psi2Sigma <- function(Psi, nu) {
  return(Psi/(nu - ncol(Psi) - 1))
}

#' @rdname RCMmisc
#' @param Sigma A numeric square positive semi-definite matrix. The expected 
#'   covariance matrix in the RCM.
#' @export
Sigma2Psi <- function(Sigma, nu) {
  return((nu - ncol(Sigma) - 1)*Sigma)
}

#' Estimate degrees of freedom 
#' 
#' Function for estimating the degrees of freedom \eqn{nu}{\nu} in the 
#' random covariance model (RCM).
#' 
#' @param Psi A numeric matrix of size \eqn{p} times \eqn{p} giving the initial
#'   estimate of \eqn{Psi}{\Psi}.
#' @param S_list A \code{list} of scatter matrices.
#' @param ns Vector of group sizes.
#' @param \dots arguments passed to the optimizer.
#' @return A list giving the \eqn{nu}{\nu} optimizing the RCM
#'   likelihood with fixed \eqn{Psi}{\Psi} and other stuff.
#' @author Anders Ellern Bilgrau
#' @note \code{rcm_get_nu} is a wrapper for \code{rcm_get_nu_optimize}.
#' @examples
#' p <- 3
#' Psi <- diag(p)
#' ns <- c(5, 5, 5)
#' true.nu <- 7
#' nus <- seq(p - 1 + sqrt(.Machine$double.eps), 40, l = 1000)
#' S_list <- createRCMData(ns = ns, psi = Psi, nu = true.nu)
#' eval.ll <- c()
#' for (i in seq_along(nus)) {
#'   eval.ll[i] <- correlateR:::rcm_loglik_arma(Psi, nus[i], S_list, ns)
#' }
#' plot(nus, eval.ll, type = "l")
#' abline(v = true.nu, col = "red", lwd = 2)
#' 
#' # Get nu
#' print(ans <- correlateR:::rcm_get_nu_optimize(Psi, S_list, ns))
#' print(ans2 <- correlateR:::rcm_get_nu_optim(Psi, S_list, ns))
#' print(ans3 <- correlateR:::rcm_get_nu_nlm(Psi, S_list, ns))
#' 
#' abline(v = ans$maximum, col = "orange", lwd = 2, lty = 2)
#' abline(v = ans2$par, col = "blue", lwd = 2, lty = 3)
#' abline(v = ans3$estimate, col = "red", lwd = 2, lty = 4)
#' 
#' \dontrun{
#' library("microbenchmark")
#' microbenchmark(correlateR:::rcm_get_nu_optimize(Psi, S_list, ns),
#'                correlateR:::rcm_get_nu_optim(Psi, S_list, ns),
#'                correlateR:::rcm_get_nu_nlm(Psi, S_list, ns))
#' }
#' @keywords internal
rcm_get_nu <- function(Psi, S_list, ns) {
  return(rcm_get_nu_optimize(Psi, S_list, ns)$maximum)
} 

#' @rdname rcm_get_nu
#' @note \code{rcm_get_nu_optimize} optimizes via \code{\link{optimize}}.
rcm_get_nu_optimize <- function(Psi, S_list, ns) {
  loglik_of_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    return(rcm_loglik_arma(Psi, nu, S_list, ns))
  }
  upper <- 1e7
  interval <- c(nrow(Psi) + 1 + sqrt(.Machine$double.eps), upper) 
  ans <- optimize(f = loglik_of_nu, interval = interval, maximum = TRUE)
  if (isTRUE(all.equal(ans$maximum, upper))) {
    stop("maximum is close to the upper edge of", upper)
  }
  return(ans)
} 

#' @rdname rcm_get_nu
#' @note \code{rcm_get_nu_optim} optimizes via \code{\link{optim}} and the L-BFGS-B
#'   method.
rcm_get_nu_optim <- function(Psi, S_list, ns, ...) {
  loglik_of_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    return(-1*rcm_loglik_arma(Psi, nu, S_list, ns))
  }
  st <- rcm_get_nu_optimize(Psi, S_list, ns)$maximum
  ans <- optim(fn = loglik_of_nu, par = st, lower = nrow(Psi) + 1, 
               method = "L-BFGS-B", hessian = TRUE, ...)
  if (ans$convergence != 0) {
    warning("optim had convergence problems; code: ", ans$convergence)
  }
  return(ans)
} 

#' @rdname rcm_get_nu
#' @note \code{rcm_get_nu_nlm} optimizes via \code{\link{nlm}}.
rcm_get_nu_nlm <- function(Psi, S_list, ns, ...) {
  loglik_of_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    return(-1*rcm_loglik_arma(Psi, nu, S_list, ns))
  }
  st <- rcm_get_nu_optimize(Psi, S_list, ns)$maximum
  ans <- nlm(f = loglik_of_nu, p = st, hessian = TRUE, ...)
  if (!(ans$code %in% c(1,2))) {
    warning("nlm had convergence problems; code: ", ans$code)
  }
  return(ans)
} 


# Compute new Psi from nu, S, ns using approximate MLE
rcm_mle_step <- function(nu, S_list, ns, ...) {
  w <- nu + ns
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) w[i]*S_list[[i]]))/sum(ns)
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
#' @param method A character giving the method to be used or abbreviation 
#'   hereof.
#' @param conf.lvl The confidence level. Default is 0.95.
#' @param eps The convergence criterion.
#' @param verbose If true, the differences in log-likelihood for each iteration
#'   is printed out.
#' @return A named list of length 3 with the elements:
#'   \item{Psi}{A matrix giving the estimate of \eqn{Psi}{Psi}.}
#'   \item{nu}{A number giving the estimate of \eqn{nu}{nu}.}
#'   \item{iterations}{A integer giving the number of iterations used.}
#'   \item{loglik}{A numeric giving the value of the log-likelihood in the last 
#'      iteration.}
#' @examples
#' nss <- c(40, 20, 30, 10)
#' print(Psii <- drop(rwishart(1)))  # Expected covariance
#' nuu <- 7
#' SS <- createRCMData(ns = nss, psi = Psii, nu = nuu)
#' 
#' fit.rcm(SS, nss, method = "EM", verbose = TRUE)
#' fit.rcm(SS, nss, method = "pool", verbose = TRUE)
#' fit.rcm(SS, nss, method = "approxMLE", verbose = FALSE)
#' @export
fit.rcm <- function(S,
                    ns,
                    Psi.init,
                    nu.init,
                    method = c("EM", "pooled", "approxMLE"),
                    conf.lvl = 0.95,
                    max.ite = 1000, 
                    eps = 1e-3,
                    verbose = FALSE) {
  method <- match.arg(method)
  p <- nrow(S[[1]])
  stopifnot(sum(ns) > p)
  
  if (missing(nu.init)) {
    nu.init <- sum(ns) + 1 + sqrt(.Machine$double.eps)
  }
  if (missing(Psi.init)) {
    Psi.init <- (nu.init - p - 1)*pool(S_list = S, ns = ns, norm_type = 1L)
  }
  
  Psi.old <- Psi.init
  nu.old <- nu.init
  ll.old <- rcm_loglik_arma(Psi = Psi.old, nu = nu.old, S_list = S, ns = ns)

  if (method == "EM" || method == "approxMLE") {
    
    updatePsi <- 
      switch(method, "EM" = rcm_em_step_arma, "approxMLE" = rcm_mle_step)

    for (i in seq_len(max.ite)) {
      Psi.new <- updatePsi(Psi = Psi.old, nu = nu.old, S_list = S, ns = ns)
      nu.new <- rcm_get_nu(Psi = Psi.new, S_list = S, ns = ns)
      ll.new <- rcm_loglik_arma(Psi = Psi.new, nu = nu.new, S_list = S, ns=ns)
      diff <- ll.new - ll.old
      if (verbose) {
        cat(sprintf("it = %g | L = %.3f | diff = %.3e\n", i, ll.new, diff))
      }
      if (diff > eps) {
        Psi.old <- Psi.new
        nu.old <- nu.new
        ll.old <- ll.new
      } else {
        if (diff < 0) {
          warning("Terminated with loglik difference smaller than 0!")
        }
        break
      }
    }

  } else if (method == "pooled") {
    
    pooled <- pool(S_list = S, ns = ns)
    
    for (i in seq_len(max.ite)) {
      nu.new <- rcm_get_nu(Psi = Psi.old, S_list = S, ns = ns)
      Psi.new <- (nu.new - p - 1)*pooled
      ll.new <- rcm_loglik_arma(Psi = Psi.old, nu = nu.new, S_list = S, ns = ns)
      diff <- ll.new - ll.old
      if (verbose) {
        cat(sprintf("it = %g | L = %.3f | diff = %.3e\n", i, ll.new, diff))
      }
      if (abs(diff) > eps) {
        Psi.old <- Psi.new
        nu.old <- nu.new
        ll.old <- ll.new
      } else {
        if (diff < 0) {
          warning("Terminated with loglik difference smaller than 0!")
        }
        break
      }
    }
      
  } else {
   
    stop("method", method, "not found.")
  
  }
  
  ans <- list("Psi" = Psi.new, "nu" = nu.new,
              "iterations" = i, "loglik" = ll.new)
  
  if (i == max.ite) {warning("max iterations (", max.ite, ") hit!")}
  return(ans)
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
  sPsi <- (nu - p - 1)*Psi
  stopifnot(ncol(x) == length(mu))
  t1 <- lgammap((nu + 1)/2, p = p)
  t2 <- -lgammap(nu/2, p = p)
  t3 <- p/2*log(pi)
  t4 <- -1/2*logdet_arma(sPsi)[1]
  
  x.center <- t(t(x) - mu)
  t5 <- -(nu + 1)/2 * log(1 + Q(x.center, solve(sPsi)))
  
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
#' @param psi The underlying \eqn{Psi} parameter which is the mean 
#'   covariance matrix. If \code{nu} is \code{Inf} 
#'   this is used as \eqn{Sigma} parameter in all Wishart distributions.
#' @param nu A numeric of length one giving the underlying \eqn{nu} parameter.
#' @param sigma A numeric matrix giving the expected covariance matrix. 
#'   Can be supplied instead of \code{psi}.
#' @return A \code{list} of matrices of the same size as \code{psi} giving
#'   observed scatter matrices from the RCM.\cr
#'   The realized covariance matrices are appended as an attribute.
#' @author Anders Ellern Bilgrau <abilgrau (at) math.aau.dk>
#' @examples
#' ns <- c(20, 10)
#' psi <- diag(3)
#' createRCMData(ns, psi, nu = 10)
#' 
#' createRCMData(ns, psi, nu = 1e20)
#' createRCMData(ns, psi, nu = Inf)
#' @export
createRCMData <- function(ns, psi, nu, sigma) {
  stopifnot(length(ns) > 0)
  k <- length(ns)
  p <- nrow(psi)
  if (nu < p + 1) {
    warning("The expected value, sigma, does not exist for n < p + 1.")
  }
  if (missing(psi) & !missing(sigma)) {
    psi <- (nu - p - 1)*sigma
  }
  if (nu == Inf) {
    sigmas <- replicate(k, psi)
  } else {
    sigmas <- rinvwishart(n = k, psi = psi, nu = nu)
  }
  S <- lapply(seq_len(k), function(i) drop(rwishart(1, sigmas[, , i], ns[i])))
  for (i in seq_along(S)) {
    dimnames(S[[i]]) <- dimnames(psi)
  }
  attributes(S)$sigmas <- sigmas
  return(S)
}


