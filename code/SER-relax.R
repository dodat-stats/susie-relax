# SER-relax: single-effect regression with relaxed LD columns.
#
# This is the main single-effect implementation from
# tex/susie-relax-report.tex. It fixes d_j = 1 and omega_j = nu0 + 3, so
# z_j | nu0 ~ N(mu_j, S_j / (nu0 + 3)).

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

make_precomp_lite <- function(Rbar) {
  J <- ncol(Rbar)
  if (max(abs(diag(Rbar) - 1)) > 1e-8) {
    stop("SER-relax fast precomputation assumes Rbar is a correlation matrix.")
  }
  Omega <- solve(Rbar)
  logdetRbar <- logdet_spd(Rbar)
  list(
    J = J,
    Omega = Omega,
    logdetRbar = logdetRbar,
    logdetS = rep(logdetRbar, J)
  )
}

add_x_precomp_lite <- function(pre, x) {
  pre$xOx <- as.numeric(crossprod(x, pre$Omega %*% x))
  pre$x_contrast_quad <- pmax(pre$xOx - x^2, 0)
  pre$Omega <- NULL
  pre
}

compute_A_B_Q_lite <- function(x, j, z_scale, z_denom, pre) {
  J <- pre$J
  quad <- pre$x_contrast_quad[j]
  Q <- (J - 1) / z_denom + z_scale^2 * quad
  B <- x[j] + z_scale * quad
  list(A = 1 + Q, B = B, Q = Q)
}

local_elbo_lite <- function(x, N, sigma0, nu0, state_j, pre, j) {
  J <- pre$J
  tau_inv2 <- nu0 + 3
  A_B_Q <- compute_A_B_Q_lite(x, j, state_j$z_scale, state_j$z_denom, pre)
  A <- A_B_Q$A
  B <- A_B_Q$B
  Q <- A_B_Q$Q

  beta_bar <- state_j$mbeta
  beta2_bar <- state_j$vbeta + state_j$mbeta^2

  eloglik <- -J / 2 * log(2 * pi) -
    0.5 * (pre$logdetRbar - J * log(N)) -
    N / 2 * (pre$xOx - 2 * beta_bar * B + beta2_bar * A)

  elog_pbeta <- -0.5 * log(2 * pi * sigma0^2) - beta2_bar / (2 * sigma0^2)

  elog_pz <- -(J - 1) / 2 * log(2 * pi) -
    0.5 * pre$logdetS[j] +
    (J - 1) / 2 * log(tau_inv2) -
    0.5 * tau_inv2 * Q

  elog_qbeta <- -0.5 * log(2 * pi * state_j$vbeta) - 0.5
  logdet_Vz <- pre$logdetS[j] - (J - 1) * log(state_j$z_denom)
  elog_qz <- -(J - 1) / 2 * log(2 * pi) -
    0.5 * logdet_Vz - (J - 1) / 2

  -log(J) + eloglik + elog_pbeta + elog_pz - elog_qbeta - elog_qz
}

update_nu0_lite <- function(states, pi, J, lower = 0.05, upper = 500) {
  Qbar <- sum(pi * vapply(states, `[[`, numeric(1), "Q"))
  tau2_hat <- Qbar / (J - 1)
  nu0_hat <- 1 / tau2_hat - 3
  min(max(nu0_hat, lower), upper)
}

ser_relax_lite_elbo <- function(local_elbo, pi) {
  pi <- as.numeric(pi)
  active <- pi > 0
  sum(pi[active] * (local_elbo[active] - log(pi[active])))
}

initialize_states_lite <- function(x, N, sigma0, nu0, pre) {
  J <- pre$J
  tau_inv2 <- nu0 + 3

  lapply(seq_len(J), function(j) {
    z_scale <- 0
    z_denom <- tau_inv2
    A_B_Q <- compute_A_B_Q_lite(x, j, z_scale, z_denom, pre)
    vbeta <- 1 / (1 / sigma0^2 + N * A_B_Q$A)
    mbeta <- vbeta * N * A_B_Q$B

    list(
      mbeta = mbeta,
      vbeta = vbeta,
      z_scale = z_scale,
      z_denom = z_denom,
      Q = A_B_Q$Q
    )
  })
}

fit_ser_relax <- function(x, Rbar, N, sigma0 = 0.2, nu0_init = 20,
                          estimate_nu0 = TRUE, max_iter = 200,
                          tol = 1e-6, verbose = FALSE,
                          nu0_bounds = c(0.05, 500),
                          nu0_update_interval = 1) {
  x <- as.numeric(x)
  Rbar <- as.matrix(Rbar)
  J <- length(x)
  stopifnot(
    nrow(Rbar) == J, ncol(Rbar) == J,
    N > 0, sigma0 > 0, nu0_init > 0,
    nu0_update_interval >= 1,
    length(nu0_bounds) == 2,
    nu0_bounds[1] > 0,
    nu0_bounds[2] > nu0_bounds[1]
  )

  pre <- add_x_precomp_lite(make_precomp_lite(Rbar), x)
  nu0 <- nu0_init
  states <- initialize_states_lite(x, N, sigma0, nu0, pre)
  pi <- rep(1 / J, J)
  elbo <- numeric()
  local <- rep(NA_real_, J)

  for (iter in seq_len(max_iter)) {
    tau_inv2 <- nu0 + 3

    for (j in seq_len(J)) {
      st <- states[[j]]

      A_B_Q <- compute_A_B_Q_lite(x, j, st$z_scale, st$z_denom, pre)
      vbeta <- 1 / (1 / sigma0^2 + N * A_B_Q$A)
      mbeta <- vbeta * N * A_B_Q$B
      beta2 <- vbeta + mbeta^2

      denom <- N * beta2 + tau_inv2
      z_scale <- N * mbeta / denom

      A_B_Q <- compute_A_B_Q_lite(x, j, z_scale, denom, pre)

      states[[j]] <- list(
        mbeta = mbeta,
        vbeta = vbeta,
        z_scale = z_scale,
        z_denom = denom,
        Q = A_B_Q$Q
      )
    }

    local <- vapply(seq_len(J), function(j) {
      local_elbo_lite(x, N, sigma0, nu0, states[[j]], pre, j)
    }, numeric(1))
    pi <- softmax(local)

    if (estimate_nu0 && iter %% nu0_update_interval == 0) {
      nu0_new <- update_nu0_lite(
        states, pi, J,
        lower = nu0_bounds[1],
        upper = nu0_bounds[2]
      )
      if (is.finite(nu0_new)) {
        nu0 <- nu0_new
      }
    }

    local <- vapply(seq_len(J), function(j) {
      local_elbo_lite(x, N, sigma0, nu0, states[[j]], pre, j)
    }, numeric(1))
    pi <- softmax(local)
    elbo[iter] <- ser_relax_lite_elbo(local, pi)

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
    tau2 = 1 / (nu0 + 3),
    local_elbo = local,
    elbo = elbo,
    states = states,
    sigma0 = sigma0,
    N = N,
    Rbar = Rbar,
    x = x
  )
}
