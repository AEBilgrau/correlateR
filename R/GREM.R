# Graphical Random Effects Model
grem_get_nu <- function(Psi, nu, S, ns, interval) {
  # Find maxima with optimize
  loglik_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    grem_loglik_nu_arma(Psi, nu, S, ns)
  }
  interval <- c(nrow(Psi) - 1 + 1e-10, 1e6)
  res <- optimize(f = loglik_nu, interval = interval, maximum = TRUE)$maximum
  return(res)
} 

# Compute new Psi from nu, S, ns using moment estimate
grem_moment_step <- function(nu, S, ns) {
  k <- length(ns)
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) S[[i]]/ns[i]))/k
  p <- nrow(Psi)
  fac <- nu - p - 1
  return(fac*Psi)
}

# Compute new Psi from nu, S, ns using approximate MLE
grem_mle_step <- function(nu, S, ns) {
  n.tot <- sum(ns)
  fac <- nu + ns
  Psi <- Reduce("+", lapply(seq_along(ns), function(i) fac[i]*S[[i]]))/n.tot
  return(Psi)
}

# Fit using the EM algorithm
#' @export
fit.grem <- function(S,
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
    ll.old  <- grem_loglik_arma(Psi.old, nu.old, S, ns)
    Psi.new <- grem_em_step_arma(Psi.old, nu.old, S, ns)
    nu.new  <- grem_get_nu(Psi.new, nu.old, S, ns, interval)
    ll.new  <- grem_loglik_arma(Psi.new, nu.new, S, ns)
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
fit.grem.MLE <- function(S, ns,
                         nu.init = nrow(S[[1]]) + 2,
                         max.ite = 1000, eps = 1e-3,
                         verbose = FALSE) {
  p <- nrow(S)
  interval <- c(p - 1 + 10*.Machine$double.eps, 1e6)
  nu.old <- nu.init
  Psi.old <- grem_mle_step(nu = nu.old, S = S, ns = ns)
  for (i in seq_len(max.ite)) {
    ll.old <- grem_loglik_arma(Psi.old, nu.old, S, ns)
    nu.new  <- grem_get_nu(Psi.old, nu.old, S, ns, interval)
    Psi.new <- grem_mle_step(nu.new, S, ns)
    ll.new  <- grem_loglik_arma(Psi.new, nu.new, S, ns)
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
fit.grem.moment <- function(S, ns,
                            nu.init = nrow(S[[1]]) + 10,
                            max.ite = 1000, eps = 1e-3,
                            verbose = FALSE) {
  p <- nrow(S[[1]])
  interval <- c(p - 1 + 10*.Machine$double.eps, 1e6)
  nu.old   <- nu.init
  Psi.old  <- grem_moment_step(nu = nu.old, S = S, ns = ns)
  for (i in seq_len(max.ite)) {
    nu.new  <- grem_get_nu(Psi = Psi.old, nu = nu.old, S = S, ns = ns,
                           interval = interval)
#     fac <- (nu.new - p - 1)/(nu.old - p - 1)
#     Psi.new <- Psi.old * fac
    Psi.new  <- grem_moment_step(nu = nu.old, S = S, ns = ns)
    if (verbose) {
      cat("ite =", i, ":", "nu.new - nu.old =", nu.new - nu.old,
          "fac =", fac, "nu =", nu.new, "\n");
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

# Density of the GREM model
#' @export
dgrem <- function(x, mu, Psi, nu, logarithm = FALSE) {
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
