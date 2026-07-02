# Collapsed-R single-effect regression.
#
# This implements the L = 1 collapsed likelihood approximation derived in
# tex/collapse-R-approach.tex. For each candidate value of nu0, the code runs
# coordinate-ascent updates for q(gamma), q_j(beta), and one of three omega
# treatments: local q_j(omega), shared q(omega), or omega = nu0 + 3. The grid
# fit selects the nu0 value with the largest optimized lower bound/objective.

logdet_spd <- function(A) {
  as.numeric(determinant(A, logarithm = TRUE)$modulus)
}

softmax <- function(x) {
  z <- x - max(x)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

log_bessel_k <- function(x, nu) {
  log(besselK(x, nu, expon.scaled = TRUE)) - x
}

small_chi <- function(lambda, chi) {
  chi <= pmax(1e-10, 1e-4 * lambda^2)
}

gig_inv_mean <- function(lambda, chi) {
  out <- numeric(length(chi))
  small <- small_chi(lambda, chi)
  out[small] <- 1 / (2 * (lambda - 1))

  if (any(!small)) {
    x <- sqrt(chi[!small])
    log_ratio <- log_bessel_k(x, lambda - 1) - log_bessel_k(x, lambda)
    vals <- exp(log_ratio - 0.5 * log(chi[!small]))
    bad <- !is.finite(vals)
    if (any(bad)) {
      vals[bad] <- 1 / (2 * (lambda - 1))
    }
    out[!small] <- vals
  }

  out
}

gig_chisq_kl <- function(lambda, chi, eta) {
  out <- numeric(length(chi))
  small <- small_chi(lambda, chi)

  if (any(!small)) {
    x <- sqrt(chi[!small])
    logK <- log_bessel_k(x, lambda)
    vals <- lgamma(lambda) +
      (lambda - 1) * log(2) -
      lambda / 2 * log(chi[!small]) -
      logK -
      chi[!small] * eta[!small] / 2
    vals[!is.finite(vals)] <- 0
    out[!small] <- vals
  }

  pmax(out, 0)
}

collapse_r_log_p0 <- function(x, Rbar, N, nu0, pre = NULL) {
  J <- length(x)
  if (is.null(pre)) {
    Omega <- solve(Rbar)
    r0 <- as.numeric(crossprod(x, Omega %*% x))
    logdetRbar <- logdet_spd(Rbar)
  } else {
    r0 <- pre$r0
    logdetRbar <- pre$logdetRbar
  }

  logdetSigma <- J * log(nu0 / (N * (nu0 + 2))) + logdetRbar
  lgamma((nu0 + J + 2) / 2) -
    lgamma((nu0 + 2) / 2) -
    J / 2 * log((nu0 + 2) * pi) -
    0.5 * logdetSigma -
    (nu0 + J + 2) / 2 * log1p(N * r0 / nu0)
}

fit_ser_collapse_r_fixed_nu0 <- function(x, Rbar, N, sigma0 = 0.2, nu0,
                                         variant = c("local", "shared", "plugin"),
                                         pi_prior = NULL, max_iter = 200,
                                         tol = 1e-6, verbose = FALSE,
                                         pre = NULL) {
  variant <- match.arg(variant)
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, sigma0 > 0, nu0 > 0,
    max_iter >= 1, tol > 0
  )

  if (is.null(pi_prior)) {
    pi_prior <- rep(1 / J, J)
  }
  pi_prior <- as.numeric(pi_prior)
  pi_prior <- pi_prior / sum(pi_prior)
  if (any(pi_prior <= 0)) {
    stop("All prior weights must be positive.")
  }

  if (is.null(pre)) {
    Omega <- solve(Rbar)
    pre <- list(
      r0 = as.numeric(crossprod(x, Omega %*% x)),
      logdetRbar = logdet_spd(Rbar)
    )
  }

  sigb <- sigma0^2
  a <- nu0 * diag(Rbar) + N * x^2
  lambda <- (nu0 + 3) / 2
  eta <- if (variant == "plugin") 1 / (nu0 + 3) else rep(1 / (nu0 + 1), J)
  alpha <- pi_prior
  elbo <- numeric()

  mu <- s2 <- m2 <- F <- rep(NA_real_, J)
  chi <- kl <- rep(NA_real_, J)
  shared_chi <- shared_kl <- NA_real_
  log_p0 <- collapse_r_log_p0(x, Rbar, N, nu0, pre)

  for (iter in seq_len(max_iter)) {
    s2 <- 1 / (1 / sigb + N * a * eta)
    mu <- s2 * N * x
    m2 <- mu^2 + s2

    if (variant == "local") {
      chi <- N * a * m2
      eta <- gig_inv_mean(lambda, chi)
      kl <- gig_chisq_kl(lambda, chi, eta)

      F <- 0.5 * (log(s2 / sigb) + 1 - m2 / sigb) +
        N * x * mu -
        0.5 * N * a * m2 * eta -
        kl

      alpha <- softmax(log(pi_prior) + F)
      omega_bound_term <- 0
    } else if (variant == "shared") {
      F <- 0.5 * (log(s2 / sigb) + 1 - m2 / sigb) +
        N * x * mu -
        0.5 * N * a * m2 * eta

      alpha <- softmax(log(pi_prior) + F)
      shared_chi <- N * sum(alpha * a * m2)
      eta <- gig_inv_mean(lambda, shared_chi)
      shared_kl <- gig_chisq_kl(lambda, shared_chi, eta)

      F <- 0.5 * (log(s2 / sigb) + 1 - m2 / sigb) +
        N * x * mu -
        0.5 * N * a * m2 * eta

      alpha <- softmax(log(pi_prior) + F)
      shared_chi <- N * sum(alpha * a * m2)
      eta <- gig_inv_mean(lambda, shared_chi)
      shared_kl <- gig_chisq_kl(lambda, shared_chi, eta)

      F <- 0.5 * (log(s2 / sigb) + 1 - m2 / sigb) +
        N * x * mu -
        0.5 * N * a * m2 * eta

      alpha <- softmax(log(pi_prior) + F)
      chi <- rep(shared_chi, J)
      kl <- rep(shared_kl, J)
      omega_bound_term <- -shared_kl
    } else {
      eta <- 1 / (nu0 + 3)
      F <- 0.5 * (log(s2 / sigb) + 1 - m2 / sigb) +
        N * x * mu -
        0.5 * N * a * m2 * eta

      alpha <- softmax(log(pi_prior) + F)
      chi <- rep(NA_real_, J)
      kl <- rep(0, J)
      omega_bound_term <- 0
    }

    active <- alpha > 0
    elbo[iter] <- log_p0 +
      omega_bound_term +
      sum(alpha[active] * (log(pi_prior[active]) - log(alpha[active]) + F[active]))

    if (verbose && (iter == 1 || iter %% 20 == 0)) {
      message(sprintf(
        "variant=%s nu0=%.4f iter=%d elbo=%.6f top=%d pi=%.4f",
        variant, nu0, iter, elbo[iter], which.max(alpha), max(alpha)
      ))
    }

    elbo_diff <- abs(elbo[iter] - elbo[iter - 1])
    if (iter > 1 && is.finite(elbo_diff) &&
        elbo_diff < tol * (1 + abs(elbo[iter - 1]))) {
      elbo <- elbo[seq_len(iter)]
      break
    }
  }

  names(alpha) <- paste0("j", seq_len(J))
  list(
    pi = alpha,
    gamma_hat = unname(which.max(alpha)),
    nu0 = nu0,
    variant = variant,
    local_elbo = F,
    elbo = elbo,
    lower_bound = tail(elbo, 1),
    states = data.frame(
      mbeta = mu,
      vbeta = s2,
      beta2 = m2,
      lambda = lambda,
      chi = chi,
      eta = eta,
      omega_kl = kl
    ),
    shared_chi = shared_chi,
    shared_omega_kl = shared_kl,
    sigma0 = sigma0,
    N = N,
    Rbar = Rbar,
    x = x
  )
}

fit_ser_collapse_r <- function(x, Rbar, N, sigma0 = 0.2,
                               nu0_grid = exp(seq(log(5), log(5000), length.out = 61)),
                               variant = c("local", "shared", "plugin"),
                               pi_prior = NULL, max_iter = 200, tol = 1e-6,
                               verbose = FALSE) {
  variant <- match.arg(variant)
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, sigma0 > 0,
    length(nu0_grid) >= 1,
    all(nu0_grid > 0)
  )

  Omega <- solve(Rbar)
  pre <- list(
    r0 = as.numeric(crossprod(x, Omega %*% x)),
    logdetRbar = logdet_spd(Rbar)
  )

  fits <- lapply(sort(unique(as.numeric(nu0_grid))), function(nu0) {
    fit_ser_collapse_r_fixed_nu0(
      x = x,
      Rbar = Rbar,
      N = N,
      sigma0 = sigma0,
      nu0 = nu0,
      variant = variant,
      pi_prior = pi_prior,
      max_iter = max_iter,
      tol = tol,
      verbose = verbose,
      pre = pre
    )
  })

  bounds <- vapply(fits, `[[`, numeric(1), "lower_bound")
  best <- which.max(bounds)
  fit <- fits[[best]]
  fit$nu0_grid <- vapply(fits, `[[`, numeric(1), "nu0")
  fit$grid_lower_bound <- bounds
  fit$all_grid_fits <- fits
  fit
}

sample_inverse_wishart <- function(df, scale) {
  W <- stats::rWishart(1, df = df, Sigma = solve(scale))[, , 1]
  solve(W)
}
