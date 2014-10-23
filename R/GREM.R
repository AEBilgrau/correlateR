# Graphical Random Effects Model

# Log multivariate gamma
grem_lgammap_old <- function(x, p = 1) {
  stopifnot(length(p) == 1)
  const <- log(pi)*p*(p - 1)/4
  terms <- sapply(seq_len(p) - 1, function(j) lgamma(x - j/2))
  dim(terms) <- c(length(x), p)
  return(const + rowSums(terms)) 
}

# Compute the log determinant easily
grem_logdet_old <- function(x, ...) {
  z <- determinant(x, logarithm = TRUE, ...)
  stopifnot(z$sign == 1)
  return(z$modulus)
}

# Log-likelihood of GREM model
grem_loglik_old <- function(Psi, nu, S, ns) {
  stopifnot(nu >= nrow(Psi) - 1)
  stopifnot(length(nu) == 1)
  
  k <- length(S)
  p <- nrow(S[[1]])
  
  cs <- (nu + ns)/2
  logdetPsi <- grem_logdet_old(Psi)[1]
  logdetPsiPlusS <- sapply(S, function(s) grem_logdet_old(Psi + s)[1])
  
  const <- sum(((ns * p)/2) * log(2))
  t1 <- k*nu/2 * logdetPsi
  t2 <- sum(grem_lgammap_old(cs, p = p))
  t3 <- -sum(cs * logdetPsiPlusS)
  t4 <- -k*grem_lgammap_old(nu/2, p = p)
  
  return(const + t1 + t2 + t3 + t4)
}


grem_loglik_nu <- function(Psi, nu, S, ns) {}
  k <- length(S)
  p <- nrow(Psi)
  logdetPsi <- correlateR:::logdet_arma(Psi)[1]
  logdetPsiplusSi <- sapply(S, function(s) correlateR:::logdet_arma(Psi + s)[1])
  t1 <- nu/2*(k*logdetPsi - sum(logdetPsiplusSi))
  t2 <- sum(correlateR:::lgammap((nu + ns)/2, p = p))
  t3 <- -k*correlateR:::lgammap(nu/2, p = p)
  return(t1 + t2 + t3)
}

grem_get_nu_new2 <- function(Psi, nu, S, ns, interval) {
  # Find maxima with optimize
  loglik_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    grem_loglik_nu(Psi, nu, S, ns)
  }
  interval <- c(nrow(Psi) - 1 + 1e-10, 1e6)
  res <- optimize(f = loglik_nu, interval = interval, maximum = TRUE)$maximum
  return(res)
} 

grem_get_nu_new3 <- function(Psi, nu, S, ns, interval) {
  # Find maxima with optimize
  loglik_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    grem_loglik_nu_arma(Psi, nu, S, ns)
  }
  interval <- c(nrow(Psi) - 1 + 1e-10, 1e6)
  res <- optimize(f = loglik_nu, interval = interval, maximum = TRUE)$maximum
  return(res)
} 

# # Derivative of the log-likelihood of GREM model
# dloglik_old <- function(Psi, nu, S, ns) {  #
#   k <- length(S)
#   p <- nrow(S[[1]])
#   t1 <-  k/2 * grem_logdet(Psi)
#   t2 <-  1/2 * sum(sapply(ns, function(ni) digammap((nu + ni)/2, p = p)))
#   t3 <- -1/2 * sum(sapply(S, function(s) grem_logdet_old(Psi + s)))
#   t4 <- -k/2 * digammap(nu/2, p = p)
#   return(t1 + t2 + t3 + t4)
# }

# Optimize loglik wrt nu
grem_get_nu_old <- function(Psi, nu, S, ns, interval) {
  # Find maxima with optimize
  loglik_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    grem_loglik_old(Psi, nu, S, ns)
  }
  interval <- c(nrow(Psi) - 1 + 1e-10, 1e6)
  res <- optimize(f = loglik_nu, interval = interval, maximum = TRUE)$maximum
  return(res)
}

grem_get_nu_new <- function(Psi, nu, S, ns, interval) {
  # Find maxima with optimize
  loglik_nu <- function(nu) { # log-likelihood as a function of nu, fixed Psi
    grem_loglik_arma(Psi, nu, S, ns)
  }
  interval <- c(nrow(Psi) - 1 + 1e-10, 1e6)
  res <- optimize(f = loglik_nu, interval = interval, maximum = TRUE)$maximum
  return(res)
} 

# get.nu2_old <- function(Psi, nu, S, ns, interval) {
#   # Find extrema as root
#   diff_loglik_nu <- function(nu) { # Derivative as a function of nu
#     dloglik(Psi, nu, S, ns)
#   }
#   res2 <- uniroot(f = diff_dloglik_nu, interval = interval)$root
#   return(res2)
# }


# Compute new Psi from Psi, nu, S, ns using the EM step
grem_em_step_old <- function(Psi, nu, S, ns) {
  k <- length(S)
  p <- nrow(S[[1]])
  co <- 1/(k*nu)
  t <- lapply(seq_along(S), function(i) (co*(ns[i] + nu))*solve(Psi + S[[i]]))
  Psi_new <- solve(Reduce("+", t))  # Sum the matrices in t
  return(Psi_new)
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
fit.grem <- function(S,
                     ns,
                     Psi.init = correlateR:::pool(S, ns),
                     nu.init = sum(ns) + 1,
                     max.ite = 1000, 
                     eps = 1e-3,
                     verbose = FALSE) {
  p <- nrow(S)
  interval <- c(p - 1 + 1e-5, 1e6)
  Psi.old <- Psi.init
  nu.old  <- nu.init
  for (i in seq_len(max.ite)) {
    ll.old  <- grem_loglik_arma(Psi.old, nu.old, S, ns)
    Psi.new <- grem_em_step_arma(Psi.old, nu.old, S, ns)
    nu.new  <- grem_get_nu_new(Psi.new, nu.old, S, ns, interval)
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



fit.grem.EM.old <- function(Psi.init,
                            nu.init,
                            S,
                            ns,
                            max.ite = 1000, eps = 1e-3,
                            verbose = FALSE) {
  p <- nrow(S)
  interval <- c(p - 1 + 1e-5, 1e6)
  Psi.old <- Psi.init
  nu.old <- nu.init
  for (i in seq_len(max.ite)) {
    ll.old <- grem_loglik_old(Psi.old, nu.old, S, ns)
    Psi.new <- grem_em_step_old(Psi.old, nu.old, S, ns)
    nu.new  <- grem_get_nu_old(Psi.new, nu.old, S, ns, interval)
    ll.new  <- grem_loglik_old(Psi.new, nu.new, S, ns)
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
fit.grem.MLE.old <- function(nu.init, S, ns,
                             max.ite = 1000, eps = 1e-3,
                             verbose = FALSE) {
  p <- nrow(S)
  interval <- c(p - 1 + 1e-5, 1e6)
  
  nu.old <- nu.init
  Psi.old <- grem_mle_step(nu = nu.old, S = S, ns = ns)
  
  for (i in seq_len(max.ite)) {
    ll.old <- grem_loglik_old(Psi.old, nu.old, S, ns)
    
    nu.new  <- grem_get_nu_old(Psi.old, nu.old, S, ns, interval)
    Psi.new <- grem_mle_step(nu.new, S, ns)
    
    ll.new  <- grem_loglik_old(Psi.new, nu.new, S, ns)
    
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
fit.grem.moment.old <- function(nu.init, S, ns,
                                max.ite = 1000, eps = 1e-3,
                                verbose = FALSE) {
  p <- nrow(S[[1]])
  interval <- c(p - 1 + 1e-5, 1e6)
  
  nu.old <- nu.init
  Psi.old <- grem_moment_step(nu = nu.old, S = S, ns = ns)
  
  for (i in seq_len(max.ite)) {
    nu.new  <- grem_get_nu_old(Psi.old, nu.old, S, ns, interval)
    fac     <- (nu.new - p - 1)/(nu.old - p - 1)
    Psi.new <- Psi.old * fac
    
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
dgrem <- function(x, mu, Psi, nu, logarithm = FALSE) {
  p <- nrow(Psi)
  if (is.null(dim(x))) {
    dim(x) <- c(1, p)
  }
  Q <- function(x, A) {
    rowSums(tcrossprod(x, A) * x)
  }
  stopifnot(ncol(x) == length(mu))
  t1 <- grem_lgammap_old((nu + 1)/2, p = p)
  t2 <- -grem_lgammap_old(nu/2, p = p)
  t3 <- p/2*log(pi)
  t4 <- -1/2*grem_logdet(Psi)
  
  x.center <- t(t(x) - mu)
  t5 <- -(nu + 1)/2 * log(1 + Q(x.center, solve(Psi)))
  
  ans <- t1 + t2 + t3 + t4 + t5
  if (!logarithm) {
    ans <- exp(ans)
  }
  attributes(ans) <- NULL
  return(ans)
}
