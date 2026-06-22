# SER-relax-prototype: single-effect regression with uncertainty in R.
#
# This file implements the coordinate-ascent variational updates derived in
# tex/susie-relax-report.tex for the original prototype single-effect model.

ar_correlation <- function(J, rho) {
  idx <- seq_len(J)
  rho ^ abs(outer(idx, idx, "-"))
}

logdet_spd <- function(A) {
  as.numeric(determinant(A, logarithm = TRUE)$modulus)
}

softmax <- function(x) {
  z <- x - max(x)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

make_precomp <- function(Rbar) {
  J <- ncol(Rbar)
  Omega <- solve(Rbar)
  out <- vector("list", J)
  for (j in seq_len(J)) {
    notj <- setdiff(seq_len(J), j)
    mu <- Rbar[notj, j]
    S <- Rbar[notj, notj, drop = FALSE] - tcrossprod(mu)
    out[[j]] <- list(
      notj = notj,
      mu = mu,
      Sinv = Omega[notj, notj, drop = FALSE],
      logdetS = logdet_spd(S)
    )
  }
  list(J = J, Omega = Omega, logdetRbar = logdet_spd(Rbar), comp = out)
}

add_x_precomp <- function(pre, x) {
  for (j in seq_len(pre$J)) {
    pc <- pre$comp[[j]]
    contrast <- x[pc$notj] - pc$mu * x[j]
    pre$comp[[j]]$x_contrast_quad <- as.numeric(crossprod(contrast, pc$Sinv %*% contrast))
    pre$comp[[j]]$Sinv <- NULL
  }
  pre$xOx <- as.numeric(crossprod(x, pre$Omega %*% x))
  pre
}

quadrature_d <- function(nu0, N, beta2_bar, A, beta_bar, B,
                         grid_size = 401, width = 10) {
  a_d <- (nu0 + 2) / 2
  b_d <- nu0 / 2

  log_kernel <- function(r) {
    -(a_d + 1) * r -
      b_d / exp(r) -
      (N / 2) * beta2_bar * A * exp(2 * r) +
      N * beta_bar * B * exp(r) +
      r
  }

  opt <- optimize(function(r) -log_kernel(r), interval = c(-25, 25))
  center <- opt$minimum
  r <- seq(center - width, center + width, length.out = grid_size)
  lr <- log_kernel(r)

  repeat {
    edge <- max(lr[1], lr[length(lr)])
    if (!is.finite(edge) || edge < max(lr) - 25 || width >= 40) {
      break
    }
    width <- width * 1.5
    r <- seq(center - width, center + width, length.out = grid_size)
    lr <- log_kernel(r)
  }

  lr_max <- max(lr)
  w_unnorm <- exp(lr - lr_max)
  delta <- r[2] - r[1]
  Z_scaled <- sum(w_unnorm) * delta
  w <- w_unnorm / sum(w_unnorm)
  logZ <- lr_max + log(Z_scaled)
  elogq_r <- sum(w * (lr - logZ))

  Ed <- sum(w * exp(r))
  Ed2 <- sum(w * exp(2 * r))
  Ed_inv <- sum(w * exp(-r))
  Elogd <- sum(w * r)

  list(
    mean = Ed,
    second = Ed2,
    inv = Ed_inv,
    log = Elogd,
    elogq = elogq_r - Elogd,
    logZ = logZ,
    mode_logd = center
  )
}

compute_A_B_Q <- function(x, j, z_scale, z_denom, pc, J) {
  quad <- pc$x_contrast_quad
  Q <- (J - 1) / z_denom + z_scale^2 * quad
  B <- x[j] + z_scale * quad
  list(A = 1 + Q, B = B, Q = Q)
}

local_elbo <- function(x, Rbar, N, sigma0, nu0, state_j, pre, j) {
  J <- pre$J
  pc <- pre$comp[[j]]
  A_B_Q <- compute_A_B_Q(x, j, state_j$z_scale, state_j$z_denom, pc, J)
  A <- A_B_Q$A
  B <- A_B_Q$B
  Q <- A_B_Q$Q

  beta_bar <- state_j$mbeta
  beta2_bar <- state_j$vbeta + state_j$mbeta^2
  d_bar <- state_j$d_mean
  d2_bar <- state_j$d_second
  omega_bar <- state_j$alpha_omega / state_j$beta_omega
  Elogomega <- digamma(state_j$alpha_omega) - log(state_j$beta_omega)

  eloglik <- -J / 2 * log(2 * pi) -
    0.5 * (pre$logdetRbar - J * log(N)) -
    N / 2 * (pre$xOx - 2 * beta_bar * d_bar * B + beta2_bar * d2_bar * A)

  elog_pbeta <- -0.5 * log(2 * pi * sigma0^2) - beta2_bar / (2 * sigma0^2)

  a_d <- (nu0 + 2) / 2
  b_d <- nu0 / 2
  elog_pd <- a_d * log(b_d) - lgamma(a_d) -
    (a_d + 1) * state_j$d_log - b_d * state_j$d_inv

  a_omega <- (nu0 + 3) / 2
  b_omega <- 0.5
  elog_pomega <- a_omega * log(b_omega) - lgamma(a_omega) +
    (a_omega - 1) * Elogomega - b_omega * omega_bar

  elog_pz <- -(J - 1) / 2 * log(2 * pi) - 0.5 * pc$logdetS +
    (J - 1) / 2 * Elogomega - 0.5 * omega_bar * Q

  elog_qbeta <- -0.5 * log(2 * pi * state_j$vbeta) - 0.5
  logdet_Vz <- pc$logdetS - (J - 1) * log(state_j$z_denom)
  elog_qz <- -(J - 1) / 2 * log(2 * pi) -
    0.5 * logdet_Vz - (J - 1) / 2
  elog_qomega <- state_j$alpha_omega * log(state_j$beta_omega) -
    lgamma(state_j$alpha_omega) +
    (state_j$alpha_omega - 1) * Elogomega -
    state_j$beta_omega * omega_bar

  -log(J) + eloglik + elog_pbeta + elog_pd + elog_pomega + elog_pz -
    elog_qbeta - state_j$d_elogq - elog_qomega - elog_qz
}

nu0_objective <- function(nu0, states, pi) {
  a_d <- (nu0 + 2) / 2
  b_d <- nu0 / 2
  a_omega <- (nu0 + 3) / 2
  b_omega <- 0.5

  vals <- vapply(seq_along(states), function(j) {
    st <- states[[j]]
    omega_bar <- st$alpha_omega / st$beta_omega
    Elogomega <- digamma(st$alpha_omega) - log(st$beta_omega)

    elog_pd <- a_d * log(b_d) - lgamma(a_d) -
      (a_d + 1) * st$d_log - b_d * st$d_inv
    elog_pomega <- a_omega * log(b_omega) - lgamma(a_omega) +
      (a_omega - 1) * Elogomega - b_omega * omega_bar
    elog_pd + elog_pomega
  }, numeric(1))

  sum(pi * vals)
}

update_nu0 <- function(nu0, states, pi, lower = 0.05, upper = 500) {
  opt <- optimize(
    f = function(log_nu) -nu0_objective(exp(log_nu), states, pi),
    interval = log(c(lower, upper))
  )
  exp(opt$minimum)
}

ser_relax_elbo <- function(local_elbo, pi) {
  pi <- as.numeric(pi)
  sum(pi * (local_elbo - log(pi)))
}

initialize_states <- function(x, Rbar, N, sigma0, nu0, pre,
                              quad_grid_size = 401, quad_width = 10) {
  J <- pre$J
  alpha_omega <- (nu0 + J + 2) / 2
  beta_omega <- 0.5
  omega_mean <- alpha_omega / beta_omega

  lapply(seq_len(J), function(j) {
    pc <- pre$comp[[j]]
    z_scale <- 0
    z_denom <- omega_mean
    A_B_Q <- compute_A_B_Q(x, j, z_scale, z_denom, pc, J)
    d_mean <- 1
    d_second <- if (nu0 > 2) nu0 / (nu0 - 2) else 2
    vbeta <- 1 / (1 / sigma0^2 + N * d_second * A_B_Q$A)
    mbeta <- vbeta * N * d_mean * A_B_Q$B
    beta2 <- vbeta + mbeta^2
    qd <- quadrature_d(
      nu0, N, beta2, A_B_Q$A, mbeta, A_B_Q$B,
      grid_size = quad_grid_size, width = quad_width
    )
    list(
      mbeta = mbeta,
      vbeta = vbeta,
      z_scale = z_scale,
      z_denom = z_denom,
      alpha_omega = alpha_omega,
      beta_omega = beta_omega,
      d_mean = qd$mean,
      d_second = qd$second,
      d_log = qd$log,
      d_inv = qd$inv,
      d_elogq = qd$elogq
    )
  })
}

fit_ser_relax_prototype <- function(x, Rbar, N, sigma0 = 0.2,
                                    nu0_init = 20,
                                    estimate_nu0 = TRUE,
                                    max_iter = 200,
                                    tol = 1e-6,
                                    verbose = FALSE,
                                    nu0_bounds = c(0.05, 500),
                                    nu0_update_interval = 1,
                                    quad_grid_size = 401,
                                    quad_width = 10) {
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, sigma0 > 0, nu0_init > 0,
    nu0_update_interval >= 1,
    quad_grid_size >= 51,
    quad_width > 0
  )

  pre <- add_x_precomp(make_precomp(Rbar), x)
  nu0 <- nu0_init
  states <- initialize_states(
    x, Rbar, N, sigma0, nu0, pre,
    quad_grid_size = quad_grid_size, quad_width = quad_width
  )
  pi <- rep(1 / J, J)
  elbo <- numeric()
  local <- rep(NA_real_, J)

  for (iter in seq_len(max_iter)) {
    for (j in seq_len(J)) {
      pc <- pre$comp[[j]]
      st <- states[[j]]

      A_B_Q <- compute_A_B_Q(x, j, st$z_scale, st$z_denom, pc, J)
      vbeta <- 1 / (1 / sigma0^2 + N * st$d_second * A_B_Q$A)
      mbeta <- vbeta * N * st$d_mean * A_B_Q$B
      beta2 <- vbeta + mbeta^2

      T1 <- mbeta * st$d_mean
      T2 <- beta2 * st$d_second
      omega_mean <- st$alpha_omega / st$beta_omega
      denom <- N * T2 + omega_mean
      z_scale <- N * T1 / denom

      A_B_Q <- compute_A_B_Q(x, j, z_scale, denom, pc, J)
      alpha_omega <- (nu0 + J + 2) / 2
      beta_omega <- (1 + A_B_Q$Q) / 2

      qd <- quadrature_d(
        nu0, N, beta2, A_B_Q$A, mbeta, A_B_Q$B,
        grid_size = quad_grid_size, width = quad_width
      )

      states[[j]] <- list(
        mbeta = mbeta,
        vbeta = vbeta,
        z_scale = z_scale,
        z_denom = denom,
        alpha_omega = alpha_omega,
        beta_omega = beta_omega,
        d_mean = qd$mean,
        d_second = qd$second,
        d_log = qd$log,
        d_inv = qd$inv,
        d_elogq = qd$elogq
      )
    }

    local <- vapply(seq_len(J), function(j) {
      local_elbo(x, Rbar, N, sigma0, nu0, states[[j]], pre, j)
    }, numeric(1))
    pi <- softmax(local)

    if (estimate_nu0 && iter %% nu0_update_interval == 0) {
      nu0_new <- update_nu0(nu0, states, pi, lower = nu0_bounds[1], upper = nu0_bounds[2])
      if (is.finite(nu0_new)) {
        nu0 <- nu0_new
      }
    }

    local <- vapply(seq_len(J), function(j) {
      local_elbo(x, Rbar, N, sigma0, nu0, states[[j]], pre, j)
    }, numeric(1))
    pi <- softmax(local)
    elbo[iter] <- ser_relax_elbo(local, pi)

    if (verbose && (iter == 1 || iter %% 10 == 0)) {
      message(sprintf(
        "iter=%d elbo=%.6f nu0=%.4f top=%d pi=%.4f",
        iter, elbo[iter], nu0, which.max(pi), max(pi)
      ))
    }

    elbo_diff <- abs(elbo[iter] - elbo[iter - 1])
    if (iter > 1 && is.finite(elbo_diff) &&
        elbo_diff < tol * (1 + abs(elbo[iter - 1]))) {
      elbo <- elbo[seq_len(iter)]
      break
    }
  }

  names(pi) <- paste0("j", seq_len(J))
  list(
    pi = pi,
    gamma_hat = which.max(pi),
    nu0 = nu0,
    local_elbo = local,
    elbo = elbo,
    states = states,
    sigma0 = sigma0,
    N = N,
    Rbar = Rbar,
    x = x
  )
}
