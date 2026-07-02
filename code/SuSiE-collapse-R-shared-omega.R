# SuSiE collapsed-R with one shared q(omega).
#
# This implements the multi-effect shared-omega approximation derived in
# tex/collapse-R-approach.tex. For each candidate nu0, it runs IBSS-style
# coordinate updates for L single effects and a single shared GIG q(omega),
# then selects the nu0 grid point with the largest optimized lower bound.

scr_logdet_spd <- function(A) {
  as.numeric(determinant(A, logarithm = TRUE)$modulus)
}

scr_softmax <- function(x) {
  z <- x - max(x)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

scr_log_bessel_k <- function(x, nu) {
  log(besselK(x, nu, expon.scaled = TRUE)) - x
}

scr_small_chi <- function(lambda, chi) {
  chi <= pmax(1e-10, 1e-4 * lambda^2)
}

scr_gig_inv_mean <- function(lambda, chi) {
  small <- scr_small_chi(lambda, chi)
  if (small) {
    return(1 / (2 * (lambda - 1)))
  }

  x <- sqrt(chi)
  val <- exp(
    scr_log_bessel_k(x, lambda - 1) -
      scr_log_bessel_k(x, lambda) -
      0.5 * log(chi)
  )
  if (is.finite(val)) val else 1 / (2 * (lambda - 1))
}

scr_gig_chisq_kl <- function(lambda, chi, eta) {
  if (scr_small_chi(lambda, chi)) {
    return(0)
  }

  x <- sqrt(chi)
  val <- lgamma(lambda) +
    (lambda - 1) * log(2) -
    lambda / 2 * log(chi) -
    scr_log_bessel_k(x, lambda) -
    chi * eta / 2

  if (is.finite(val)) max(val, 0) else 0
}

scr_log_p0 <- function(x, Rbar, N, nu0, pre = NULL) {
  J <- length(x)
  if (is.null(pre)) {
    Omega <- solve(Rbar)
    r0 <- as.numeric(crossprod(x, Omega %*% x))
    logdetRbar <- scr_logdet_spd(Rbar)
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

scr_initialize_effect <- function(J, sigma2, pi_prior) {
  alpha <- pi_prior
  mu <- rep(0, J)
  s2 <- rep(sigma2, J)
  m2 <- mu^2 + s2
  m <- alpha * mu

  list(
    alpha = alpha,
    mu = mu,
    s2 = s2,
    m2 = m2,
    m = m,
    score = rep(NA_real_, J)
  )
}

scr_update_effect <- function(effect, m_minus, x, A, diagA, N, sigma2,
                              eta, pi_prior) {
  cvec <- as.numeric(A %*% m_minus)
  natural_mean <- x - eta * cvec
  s2 <- 1 / (1 / sigma2 + N * eta * diagA)
  mu <- s2 * N * natural_mean
  m2 <- mu^2 + s2

  score <- 0.5 * (log(s2 / sigma2) + 1 - m2 / sigma2) +
    N * mu * natural_mean -
    0.5 * N * eta * diagA * m2
  alpha <- scr_softmax(log(pi_prior) + score)
  names(alpha) <- paste0("j", seq_along(alpha))

  list(
    alpha = alpha,
    mu = mu,
    s2 = s2,
    m2 = m2,
    m = as.numeric(alpha) * mu,
    score = score
  )
}

scr_effect_second_moment <- function(effects, A) {
  m_total <- Reduce(`+`, lapply(effects, `[[`, "m"))
  sq <- as.numeric(crossprod(m_total, A %*% m_total))
  diagA <- diag(A)

  for (eff in effects) {
    sq <- sq +
      sum(as.numeric(eff$alpha) * eff$m2 * diagA) -
      as.numeric(crossprod(eff$m, A %*% eff$m))
  }

  sq
}

scr_lower_bound <- function(x, Rbar, N, nu0, sigma2, effects, A, eta,
                            omega_kl, pi_prior, pre) {
  m_total <- Reduce(`+`, lapply(effects, `[[`, "m"))
  S_Q <- scr_effect_second_moment(effects, A)
  out <- scr_log_p0(x, Rbar, N, nu0, pre) +
    N * sum(x * m_total) -
    0.5 * N * eta * S_Q -
    omega_kl

  for (ell in seq_along(effects)) {
    eff <- effects[[ell]]
    active <- as.numeric(eff$alpha) > 0
    normal_term <- 0.5 * (
      log(eff$s2 / sigma2[ell]) + 1 - eff$m2 / sigma2[ell]
    )
    out <- out + sum(as.numeric(eff$alpha)[active] * (
      log(pi_prior[active]) - log(as.numeric(eff$alpha)[active]) +
        normal_term[active]
    ))
  }

  out
}

fit_susie_collapse_r_shared_fixed_nu0 <- function(
    x, Rbar, N, L = 2, sigma2 = 0.2^2, nu0,
    pi_prior = NULL, estimate_sigma2 = FALSE,
    max_iter = 100, tol = 1e-6, verbose = FALSE, pre = NULL) {
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  if (length(sigma2) == 1) {
    sigma2 <- rep(sigma2, L)
  }
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, L >= 1, length(sigma2) == L,
    all(sigma2 > 0), nu0 > 0,
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
      logdetRbar = scr_logdet_spd(Rbar)
    )
  }

  A <- nu0 * Rbar + N * tcrossprod(x)
  diagA <- diag(A)
  lambda <- (nu0 + 3) / 2
  eta <- 1 / (nu0 + 1)
  chi <- NA_real_
  omega_kl <- NA_real_
  effects <- lapply(seq_len(L), function(ell) {
    scr_initialize_effect(J, sigma2[ell], pi_prior)
  })
  elbo <- numeric()

  for (iter in seq_len(max_iter)) {
    m_sum <- Reduce(`+`, lapply(effects, `[[`, "m"))

    for (ell in seq_len(L)) {
      m_minus <- m_sum - effects[[ell]]$m
      old_m <- effects[[ell]]$m
      effects[[ell]] <- scr_update_effect(
        effects[[ell]], m_minus, x, A, diagA, N, sigma2[ell], eta, pi_prior
      )
      m_sum <- m_sum - old_m + effects[[ell]]$m
    }

    if (estimate_sigma2) {
      sigma2 <- vapply(effects, function(eff) {
        max(sum(as.numeric(eff$alpha) * eff$m2), 1e-12)
      }, numeric(1))
    }

    S_Q <- scr_effect_second_moment(effects, A)
    chi <- N * S_Q
    eta <- scr_gig_inv_mean(lambda, chi)
    omega_kl <- scr_gig_chisq_kl(lambda, chi, eta)
    elbo[iter] <- scr_lower_bound(
      x, Rbar, N, nu0, sigma2, effects, A, eta, omega_kl, pi_prior, pre
    )

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      tops <- vapply(effects, function(eff) which.max(eff$alpha), integer(1))
      pips <- vapply(effects, function(eff) max(eff$alpha), numeric(1))
      message(sprintf(
        "nu0=%.4f iter=%d elbo=%.6f eta=%.4g top=(%s) pi=(%s)",
        nu0, iter, elbo[iter], eta,
        paste(tops, collapse = ","),
        paste(sprintf("%.3f", pips), collapse = ",")
      ))
    }

    elbo_diff <- abs(elbo[iter] - elbo[iter - 1])
    if (iter > 1 && is.finite(elbo_diff) &&
        elbo_diff < tol * (1 + abs(elbo[iter - 1]))) {
      elbo <- elbo[seq_len(iter)]
      break
    }
  }

  alpha <- do.call(rbind, lapply(effects, function(eff) as.numeric(eff$alpha)))
  rownames(alpha) <- paste0("effect", seq_len(L))
  colnames(alpha) <- paste0("j", seq_len(J))
  pip <- 1 - apply(1 - alpha, 2, prod)

  list(
    alpha = alpha,
    pip = pip,
    gamma_hat = apply(alpha, 1, which.max),
    nu0 = nu0,
    eta = eta,
    chi = chi,
    omega_kl = omega_kl,
    lower_bound = tail(elbo, 1),
    elbo = elbo,
    effects = effects,
    sigma2 = sigma2,
    N = N,
    Rbar = Rbar,
    x = x
  )
}

fit_susie_collapse_r_shared <- function(
    x, Rbar, N, L = 2, sigma2 = 0.2^2,
    nu0_grid = exp(seq(log(5), log(5000), length.out = 61)),
    pi_prior = NULL, estimate_sigma2 = FALSE,
    max_iter = 100, tol = 1e-6, verbose = FALSE) {
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, L >= 1,
    length(nu0_grid) >= 1,
    all(nu0_grid > 0)
  )

  Omega <- solve(Rbar)
  pre <- list(
    r0 = as.numeric(crossprod(x, Omega %*% x)),
    logdetRbar = scr_logdet_spd(Rbar)
  )

  fits <- lapply(sort(unique(as.numeric(nu0_grid))), function(nu0) {
    fit_susie_collapse_r_shared_fixed_nu0(
      x = x,
      Rbar = Rbar,
      N = N,
      L = L,
      sigma2 = sigma2,
      nu0 = nu0,
      pi_prior = pi_prior,
      estimate_sigma2 = estimate_sigma2,
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

scr_sample_inverse_wishart <- function(df, scale) {
  W <- stats::rWishart(1, df = df, Sigma = solve(scale))[, , 1]
  solve(W)
}
