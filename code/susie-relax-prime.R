# SuSiE-relax-prime: multi-effect regression with latent relaxed
# contribution vectors.
#
# This implements tex/susie-relax-prime.tex. Conditional on each active
# coordinate, beta has a scalar Gaussian posterior and the relaxed component
# eta has a Gaussian posterior whose covariance is proportional to the Schur
# complement S_j. The implementation uses the identity
#   (r_-j - mu_j r_j)' S_j^{-1} (r_-j - mu_j r_j) = r' Omega r - r_j^2
# and never forms the S_j matrices.

srp_logdet_spd <- function(A) {
  as.numeric(determinant(A, logarithm = TRUE)$modulus)
}

srp_softmax <- function(x) {
  z <- x - max(x)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

make_susie_relax_prime_precomp <- function(Rbar) {
  J <- ncol(Rbar)
  if (max(abs(diag(Rbar) - 1)) > 1e-8) {
    stop("SuSiE-relax-prime fast precomputation assumes Rbar is a correlation matrix.")
  }
  list(
    J = J,
    Omega = solve(Rbar),
    logdetRbar = srp_logdet_spd(Rbar)
  )
}

susie_relax_prime_residual_stats <- function(r, pre) {
  rOr <- as.numeric(crossprod(r, pre$Omega %*% r))
  list(
    rOr = rOr,
    contrast_quad = pmax(rOr - r^2, 0)
  )
}

susie_relax_prime_effect_mean <- function(Rbar, r, alpha, m, rho) {
  beta_mean <- as.numeric(Rbar %*% (alpha * m))
  eta_mean <- rho * r - rho * as.numeric(Rbar %*% (alpha * r))
  beta_mean + eta_mean
}

initialize_susie_relax_prime_effect <- function(J, sigma2, tau2) {
  alpha <- rep(1 / J, J)
  names(alpha) <- paste0("j", seq_len(J))
  list(
    alpha = alpha,
    m = rep(0, J),
    v = sigma2,
    kappa = 0,
    rho = 0,
    c = sigma2 * tau2,
    B = rep(sigma2, J),
    D = rep((J - 1) * sigma2 * tau2, J),
    theta = rep(0, J),
    log_weights = rep(0, J)
  )
}

update_susie_relax_prime_effect <- function(r, Rbar, N, sigma2, tau2, pre) {
  J <- pre$J
  r_stats <- susie_relax_prime_residual_stats(r, pre)

  v <- 1 / (N + 1 / sigma2)
  kappa <- N * v
  rho <- N * sigma2 * tau2 / (1 + N * sigma2 * tau2)
  c <- sigma2 * tau2 / (1 + N * sigma2 * tau2)

  m <- kappa * r
  log_weights <- N * (kappa - rho) * r^2 / 2
  alpha <- srp_softmax(log_weights)
  names(alpha) <- paste0("j", seq_len(J))

  B <- m^2 + v
  D <- rho^2 * r_stats$contrast_quad + (J - 1) * c
  theta <- susie_relax_prime_effect_mean(Rbar, r, alpha, m, rho)

  list(
    alpha = alpha,
    m = m,
    v = v,
    kappa = kappa,
    rho = rho,
    c = c,
    B = B,
    D = D,
    theta = theta,
    log_weights = log_weights
  )
}

update_susie_relax_prime_effects_sweep <- function(effects, x, Rbar, N,
                                                   sigma2, tau2, pre) {
  L <- length(effects)
  theta_sum <- Reduce(`+`, lapply(effects, `[[`, "theta"))
  for (ell in seq_len(L)) {
    r <- x - theta_sum + effects[[ell]]$theta
    old_theta <- effects[[ell]]$theta
    effects[[ell]] <- update_susie_relax_prime_effect(
      r, Rbar, N, sigma2[ell], tau2, pre
    )
    theta_sum <- theta_sum - old_theta + effects[[ell]]$theta
  }
  effects
}

update_susie_relax_prime_sigma2 <- function(effects, tau2,
                                            lower = 1e-12, upper = Inf) {
  vapply(effects, function(eff) {
    Bbar <- sum(eff$alpha * eff$B)
    Dbar <- sum(eff$alpha * eff$D)
    sigma2_hat <- (Bbar + Dbar / tau2) / length(eff$alpha)
    min(max(sigma2_hat, lower), upper)
  }, numeric(1))
}

update_susie_relax_prime_nu0 <- function(effects, sigma2, J,
                                         lower = 0.05, upper = 500) {
  D_scaled <- 0
  for (ell in seq_along(effects)) {
    D_scaled <- D_scaled + sum(effects[[ell]]$alpha * effects[[ell]]$D) /
      sigma2[ell]
  }
  tau2_hat <- D_scaled / (length(effects) * (J - 1))
  nu0_hat <- 1 / tau2_hat - 3
  min(max(nu0_hat, lower), upper)
}

susie_relax_prime_elbo <- function(x, Rbar, N, sigma2, tau2, effects, pre) {
  J <- pre$J
  theta_sum <- Reduce(`+`, lapply(effects, `[[`, "theta"))
  residual <- x - theta_sum
  residual_quad <- as.numeric(crossprod(residual, pre$Omega %*% residual))

  variance_correction <- 0
  kl_sum <- 0
  for (ell in seq_along(effects)) {
    eff <- effects[[ell]]
    theta_quad <- as.numeric(crossprod(eff$theta, pre$Omega %*% eff$theta))
    expected_theta_quad <- sum(eff$alpha * (eff$B + eff$D))
    effect_variance <- expected_theta_quad - theta_quad
    if (effect_variance < 0 && effect_variance > -sqrt(.Machine$double.eps)) {
      effect_variance <- 0
    }
    variance_correction <- variance_correction + effect_variance

    active <- eff$alpha > 0
    kl_beta <- 0.5 * (
      log(sigma2[ell] / eff$v) + eff$B / sigma2[ell] - 1
    )
    kl_eta <- 0.5 * (
      (J - 1) * log(sigma2[ell] * tau2 / eff$c) +
        eff$D / (sigma2[ell] * tau2) - (J - 1)
    )
    kl_sum <- kl_sum + sum(eff$alpha[active] * (
      log(J * eff$alpha[active]) + kl_beta[active] + kl_eta[active]
    ))
  }

  Cx <- -J / 2 * log(2 * pi) -
    0.5 * (pre$logdetRbar - J * log(N))
  Cx - N / 2 * (residual_quad + variance_correction) - kl_sum
}

fit_susie_relax_prime <- function(x, Rbar, N, L = 2, sigma2 = 0.3^2,
                                  nu0_init = 1000,
                                  estimate_nu0 = TRUE,
                                  estimate_sigma2 = FALSE,
                                  warmup_iter = 5,
                                  max_iter = 100, tol = 1e-6,
                                  verbose = FALSE,
                                  nu0_bounds = c(1, 2000),
                                  sigma2_bounds = c(1e-12, Inf),
                                  elbo_update_interval = 1,
                                  nu0_update_interval = 1,
                                  sigma2_update_interval = 1,
                                  nu0_tol = Inf) {
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  if (length(sigma2) == 1) {
    sigma2 <- rep(sigma2, L)
  }
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    length(sigma2) == L,
    N > 0, L >= 1,
    all(sigma2 > 0),
    nu0_init > 0,
    length(nu0_bounds) == 2,
    nu0_bounds[1] > 0,
    nu0_bounds[2] > nu0_bounds[1],
    length(sigma2_bounds) == 2,
    sigma2_bounds[1] > 0,
    sigma2_bounds[2] >= sigma2_bounds[1],
    elbo_update_interval >= 1,
    nu0_update_interval >= 1,
    sigma2_update_interval >= 1
  )

  pre <- make_susie_relax_prime_precomp(Rbar)
  nu0 <- nu0_init
  tau2 <- 1 / (nu0 + 3)
  effects <- lapply(seq_len(L), function(ell) {
    initialize_susie_relax_prime_effect(J, sigma2[ell], tau2)
  })

  elbo <- numeric()
  nu0_trace <- numeric()
  sigma2_trace <- matrix(NA_real_, nrow = max_iter, ncol = L)
  colnames(sigma2_trace) <- paste0("effect", seq_len(L))

  for (iter in seq_len(max_iter)) {
    nu0_prev <- nu0
    effects <- update_susie_relax_prime_effects_sweep(
      effects, x, Rbar, N, sigma2, tau2, pre
    )

    can_update_sigma2 <- estimate_sigma2 && iter > warmup_iter &&
      iter %% sigma2_update_interval == 0
    if (can_update_sigma2) {
      sigma2 <- update_susie_relax_prime_sigma2(
        effects, tau2,
        lower = sigma2_bounds[1],
        upper = sigma2_bounds[2]
      )
    }

    can_update_nu0 <- estimate_nu0 && iter > warmup_iter &&
      iter %% nu0_update_interval == 0
    if (can_update_nu0) {
      nu0_new <- update_susie_relax_prime_nu0(
        effects, sigma2, J,
        lower = nu0_bounds[1],
        upper = nu0_bounds[2]
      )
      if (is.finite(nu0_new)) {
        nu0 <- nu0_new
        tau2 <- 1 / (nu0 + 3)
      }
    }

    should_compute_elbo <- iter == 1 || iter %% elbo_update_interval == 0 ||
      iter == max_iter
    if (should_compute_elbo) {
      elbo[iter] <- susie_relax_prime_elbo(
        x, Rbar, N, sigma2, tau2, effects, pre
      )
    } else {
      elbo[iter] <- NA_real_
    }
    nu0_trace[iter] <- nu0
    sigma2_trace[iter, ] <- sigma2

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      tops <- vapply(effects, function(eff) which.max(eff$alpha), integer(1))
      pips <- vapply(effects, function(eff) max(eff$alpha), numeric(1))
      elbo_text <- if (is.na(elbo[iter])) "NA" else sprintf("%.6f", elbo[iter])
      message(sprintf(
        "iter=%d elbo=%s nu0=%.3f top=(%s) pip=(%s)",
        iter, elbo_text, nu0,
        paste(tops, collapse = ","),
        paste(sprintf("%.3f", pips), collapse = ",")
      ))
    }

    previous_elbo_iter <- if (iter > 1) {
      tail(which(!is.na(elbo[seq_len(iter - 1)])), 1)
    } else {
      integer(0)
    }
    if (should_compute_elbo && length(previous_elbo_iter) == 1) {
      elbo_diff <- abs(elbo[iter] - elbo[previous_elbo_iter])
      nu0_diff <- abs(log(nu0) - log(nu0_prev))
      elbo_converged <- is.finite(elbo_diff) &&
        elbo_diff < tol * (1 + abs(elbo[previous_elbo_iter]))
      nu0_converged <- is.finite(nu0_diff) && nu0_diff < nu0_tol
      if (elbo_converged && nu0_converged) {
        elbo <- elbo[seq_len(iter)]
        nu0_trace <- nu0_trace[seq_len(iter)]
        sigma2_trace <- sigma2_trace[seq_len(iter), , drop = FALSE]
        break
      }
    }
  }

  alpha <- do.call(rbind, lapply(effects, `[[`, "alpha"))
  rownames(alpha) <- paste0("effect", seq_len(L))
  colnames(alpha) <- paste0("j", seq_len(J))
  pip <- 1 - apply(1 - alpha, 2, prod)

  list(
    alpha = alpha,
    pip = pip,
    gamma_hat = apply(alpha, 1, which.max),
    nu0 = nu0,
    tau2 = tau2,
    nu0_trace = nu0_trace,
    effects = effects,
    elbo = elbo,
    sigma2 = sigma2,
    sigma2_trace = sigma2_trace,
    N = N,
    Rbar = Rbar,
    x = x,
    warmup_iter = warmup_iter,
    elbo_update_interval = elbo_update_interval,
    nu0_tol = nu0_tol
  )
}
