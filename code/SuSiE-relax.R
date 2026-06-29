# SuSiE-relax: multi-effect regression with relaxed LD columns.
#
# This is the main multi-effect implementation from
# tex/susie-relax-report.tex, using
# d_{ell j} = 1 and omega_{ell j} = nu0 + 3. The implementation uses
# the correlation-matrix identities
#   (r_-j - mu_j r_j)' S_j^{-1} (r_-j - mu_j r_j) = r' Omega r - r_j^2
# and log |S_j| = log |Rbar|.

logdet_spd <- function(A) {
  as.numeric(determinant(A, logarithm = TRUE)$modulus)
}

softmax <- function(x) {
  z <- x - max(x)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

make_susie_relax_precomp <- function(Rbar) {
  J <- ncol(Rbar)
  if (max(abs(diag(Rbar) - 1)) > 1e-8) {
    stop("SuSiE-relax fast precomputation assumes Rbar is a correlation matrix.")
  }
  Omega <- solve(Rbar)
  list(
    J = J,
    Omega = Omega,
    logdetRbar = logdet_spd(Rbar)
  )
}

susie_relax_residual_stats <- function(r, pre) {
  rOr <- as.numeric(crossprod(r, pre$Omega %*% r))
  list(
    rOr = rOr,
    quad = pmax(rOr - r^2, 0)
  )
}

susie_relax_coordinate_stats <- function(r, z_scale, z_denom, pre, r_stats) {
  Q <- (pre$J - 1) / z_denom + z_scale^2 * r_stats$quad
  B <- r + z_scale * r_stats$quad
  list(A = 1 + Q, B = B, Q = Q)
}

susie_relax_effect_mean <- function(Rbar, r, pi, mbeta, z_scale) {
  weights <- pi * mbeta
  col_weights <- weights * (1 - z_scale * r)
  residual_weight <- sum(weights * z_scale)
  as.numeric(Rbar %*% col_weights + residual_weight * r)
}

susie_relax_local_elbo <- function(r, N, sigma2, nu0, state, pre, r_stats) {
  J <- pre$J
  tau_inv2 <- nu0 + 3
  beta2 <- state$vbeta + state$mbeta^2

  eloglik <- -J / 2 * log(2 * pi) -
    0.5 * (pre$logdetRbar - J * log(N)) -
    N / 2 * (r_stats$rOr - 2 * state$mbeta * state$B + beta2 * state$A)

  elog_pbeta <- -0.5 * log(2 * pi * sigma2) - beta2 / (2 * sigma2)

  elog_pz <- -(J - 1) / 2 * log(2 * pi) -
    0.5 * pre$logdetRbar +
    (J - 1) / 2 * log(tau_inv2) -
    0.5 * tau_inv2 * state$Q

  elog_qbeta <- -0.5 * log(2 * pi * state$vbeta) - 0.5
  logdet_Vz <- pre$logdetRbar - (J - 1) * log(state$z_denom)
  elog_qz <- -(J - 1) / 2 * log(2 * pi) -
    0.5 * logdet_Vz - (J - 1) / 2

  log(1 / J) + eloglik + elog_pbeta + elog_pz - elog_qbeta - elog_qz
}

susie_relax_prior_entropy_terms <- function(state, sigma2, nu0, pre) {
  J <- pre$J
  tau_inv2 <- nu0 + 3
  beta2 <- state$vbeta + state$mbeta^2

  elog_pbeta <- -0.5 * log(2 * pi * sigma2) - beta2 / (2 * sigma2)
  elog_pz <- -(J - 1) / 2 * log(2 * pi) -
    0.5 * pre$logdetRbar +
    (J - 1) / 2 * log(tau_inv2) -
    0.5 * tau_inv2 * state$Q

  elog_qbeta <- -0.5 * log(2 * pi * state$vbeta) - 0.5
  logdet_Vz <- pre$logdetRbar - (J - 1) * log(state$z_denom)
  elog_qz <- -(J - 1) / 2 * log(2 * pi) -
    0.5 * logdet_Vz - (J - 1) / 2

  log(1 / J) + elog_pbeta + elog_pz - elog_qbeta - elog_qz
}

initialize_susie_relax_effect <- function(x, Rbar, N, sigma2, nu0, pre) {
  J <- pre$J
  tau_inv2 <- nu0 + 3
  r_stats <- susie_relax_residual_stats(x, pre)
  z_scale <- rep(0, J)
  z_denom <- rep(tau_inv2, J)
  stats <- susie_relax_coordinate_stats(x, z_scale, z_denom, pre, r_stats)
  vbeta <- 1 / (1 / sigma2 + N * stats$A)
  mbeta <- vbeta * N * stats$B

  states <- lapply(seq_len(J), function(j) {
    list(
      mbeta = mbeta[j],
      vbeta = vbeta[j],
      z_scale = z_scale[j],
      z_denom = z_denom[j],
      A = stats$A[j],
      B = stats$B[j],
      Q = stats$Q[j]
    )
  })

  pi <- rep(1 / J, J)
  theta <- susie_relax_effect_mean(Rbar, x, pi, mbeta, z_scale)
  names(pi) <- paste0("j", seq_len(J))
  list(
    pi = pi,
    mbeta = mbeta,
    vbeta = vbeta,
    z_scale = z_scale,
    z_denom = z_denom,
    A = stats$A,
    B = stats$B,
    Q = stats$Q,
    states = states,
    theta = theta,
    local_elbo = rep(NA_real_, J)
  )
}

update_susie_relax_effect <- function(effect, r, Rbar, N, sigma2, nu0, pre) {
  J <- pre$J
  tau_inv2 <- nu0 + 3
  r_stats <- susie_relax_residual_stats(r, pre)

  stats <- susie_relax_coordinate_stats(r, effect$z_scale, effect$z_denom, pre, r_stats)
  vbeta <- 1 / (1 / sigma2 + N * stats$A)
  mbeta <- vbeta * N * stats$B
  beta2 <- vbeta + mbeta^2

  z_denom <- N * beta2 + tau_inv2
  z_scale <- N * mbeta / z_denom
  stats <- susie_relax_coordinate_stats(r, z_scale, z_denom, pre, r_stats)

  states <- lapply(seq_len(J), function(j) {
    list(
      mbeta = mbeta[j],
      vbeta = vbeta[j],
      z_scale = z_scale[j],
      z_denom = z_denom[j],
      A = stats$A[j],
      B = stats$B[j],
      Q = stats$Q[j]
    )
  })

  local <- vapply(seq_len(J), function(j) {
    susie_relax_local_elbo(r, N, sigma2, nu0, states[[j]], pre, r_stats)
  }, numeric(1))

  pi_new <- softmax(local)

  theta <- susie_relax_effect_mean(Rbar, r, pi_new, mbeta, z_scale)
  names(pi_new) <- paste0("j", seq_len(J))
  list(
    pi = pi_new,
    mbeta = mbeta,
    vbeta = vbeta,
    z_scale = z_scale,
    z_denom = z_denom,
    A = stats$A,
    B = stats$B,
    Q = stats$Q,
    states = states,
    theta = theta,
    local_elbo = local
  )
}

update_susie_relax_effects_sweep <- function(effects, x, Rbar, N, sigma2,
                                            nu0, pre) {
  L <- length(effects)
  theta_sum <- Reduce(`+`, lapply(effects, `[[`, "theta"))
  for (ell in seq_len(L)) {
    r <- x - theta_sum + effects[[ell]]$theta
    old_theta <- effects[[ell]]$theta
    effects[[ell]] <- update_susie_relax_effect(
      effects[[ell]], r, Rbar, N, sigma2[ell], nu0, pre
    )
    theta_sum <- theta_sum - old_theta + effects[[ell]]$theta
  }
  effects
}

update_susie_relax_nu0 <- function(effects, J, lower = 0.05, upper = 500) {
  total_Q <- 0
  for (eff in effects) {
    total_Q <- total_Q + sum(eff$pi * eff$Q)
  }
  tau2_hat <- total_Q / (length(effects) * (J - 1))
  nu0_hat <- 1 / tau2_hat - 3
  min(max(nu0_hat, lower), upper)
}

update_susie_relax_sigma2 <- function(effects, lower = 1e-12, upper = Inf) {
  vapply(effects, function(eff) {
    beta2 <- eff$vbeta + eff$mbeta^2
    sigma2_hat <- sum(eff$pi * beta2)
    min(max(sigma2_hat, lower), upper)
  }, numeric(1))
}

categorical_entropy <- function(pi) {
  pi <- as.numeric(pi)
  active <- pi > 0
  -sum(pi[active] * log(pi[active]))
}

susie_relax_elbo <- function(x, Rbar, N, sigma2, nu0, effects, pre) {
  J <- pre$J
  L <- length(effects)
  theta_sum <- Reduce(`+`, lapply(effects, `[[`, "theta"))
  xOx <- as.numeric(crossprod(x, pre$Omega %*% x))
  theta_x <- as.numeric(crossprod(theta_sum, pre$Omega %*% x))

  self_second <- 0
  prior_entropy <- 0
  for (ell in seq_len(L)) {
    eff <- effects[[ell]]
    beta2 <- eff$vbeta + eff$mbeta^2
    self_second <- self_second + sum(eff$pi * beta2 * eff$A)

    terms <- vapply(seq_len(J), function(j) {
      susie_relax_prior_entropy_terms(eff$states[[j]], sigma2[ell], nu0, pre)
    }, numeric(1))
    prior_entropy <- prior_entropy + sum(eff$pi * terms) + categorical_entropy(eff$pi)
  }

  cross <- 0
  if (L > 1) {
    for (ell in seq_len(L - 1)) {
      for (k in (ell + 1):L) {
        cross <- cross +
          2 * as.numeric(crossprod(effects[[ell]]$theta, pre$Omega %*% effects[[k]]$theta))
      }
    }
  }

  eloglik <- -J / 2 * log(2 * pi) -
    0.5 * (pre$logdetRbar - J * log(N)) -
    N / 2 * (xOx - 2 * theta_x + self_second + cross)

  eloglik + prior_entropy
}

fit_susie_relax <- function(x, Rbar, N, L = 2, sigma2 = 0.3^2,
                            nu0_init = 1000, estimate_nu0 = TRUE,
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

  pre <- make_susie_relax_precomp(Rbar)
  nu0 <- nu0_init
  effects <- vector("list", L)
  for (ell in seq_len(L)) {
    effects[[ell]] <- initialize_susie_relax_effect(x, Rbar, N, sigma2[ell], nu0, pre)
  }

  elbo <- numeric()
  nu0_trace <- numeric()
  sigma2_trace <- matrix(NA_real_, nrow = max_iter, ncol = L)
  colnames(sigma2_trace) <- paste0("effect", seq_len(L))
  for (iter in seq_len(max_iter)) {
    nu0_prev <- nu0
    effects <- update_susie_relax_effects_sweep(
      effects, x, Rbar, N, sigma2, nu0, pre
    )

    can_update_sigma2 <- estimate_sigma2 && iter > warmup_iter &&
      iter %% sigma2_update_interval == 0
    if (can_update_sigma2) {
      sigma2 <- update_susie_relax_sigma2(
        effects,
        lower = sigma2_bounds[1],
        upper = sigma2_bounds[2]
      )
    }

    can_update_nu0 <- estimate_nu0 && iter > warmup_iter &&
      iter %% nu0_update_interval == 0
    if (can_update_nu0) {
      nu0_new <- update_susie_relax_nu0(
        effects, J,
        lower = nu0_bounds[1],
        upper = nu0_bounds[2]
      )
      if (is.finite(nu0_new)) {
        nu0 <- nu0_new
      }
    }

    should_compute_elbo <- iter == 1 || iter %% elbo_update_interval == 0 ||
      iter == max_iter
    if (should_compute_elbo) {
      elbo[iter] <- susie_relax_elbo(x, Rbar, N, sigma2, nu0, effects, pre)
    } else {
      elbo[iter] <- NA_real_
    }
    nu0_trace[iter] <- nu0
    sigma2_trace[iter, ] <- sigma2

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      tops <- vapply(effects, function(eff) which.max(eff$pi), integer(1))
      pips <- vapply(effects, function(eff) max(eff$pi), numeric(1))
      elbo_text <- if (is.na(elbo[iter])) "NA" else sprintf("%.6f", elbo[iter])
      message(sprintf(
        "iter=%d elbo=%s nu0=%.3f sigma2=(%s) top=(%s) pip=(%s)",
        iter, elbo_text, nu0,
        paste(sprintf("%.4g", sigma2), collapse = ","),
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

  alpha <- do.call(rbind, lapply(effects, `[[`, "pi"))
  rownames(alpha) <- paste0("effect", seq_len(L))
  colnames(alpha) <- paste0("j", seq_len(J))
  pip <- 1 - apply(1 - alpha, 2, prod)

  list(
    alpha = alpha,
    pip = pip,
    gamma_hat = apply(alpha, 1, which.max),
    nu0 = nu0,
    tau2 = 1 / (nu0 + 3),
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
