# Compare SuSiE-relax and collapsed-R shared-omega SuSiE across seeds in the
# different-ancestry setting. Both methods learn lambda by grid search and learn
# nu0 inside each lambda fit.
#
# Example:
#   source("analysis/compare-collapse-vs-relax-diff-ancestry.R")
#   cmp <- run_diff_ancestry_comparison(seeds = 1:25)
#   cmp$selected_L_distribution

library(susieR)

source("code/SuSiE-relax.R")
source("code/SuSiE-collapse-R-shared-omega.R")

make_rbar <- function(R0, lambda = 0.01) {
  Rbar <- (1 - lambda) * R0 + lambda * diag(1, ncol(R0))
  cov2cor(Rbar)
}

make_cs_from_alpha <- function(alpha, R, coverage = 0.95,
                               min_abs_corr = 0.5,
                               min_top_alpha = 1e-4) {
  alpha <- as.matrix(alpha)
  cs <- list()
  purity <- numeric()

  for (ell in seq_len(nrow(alpha))) {
    a <- as.numeric(alpha[ell, ])
    if (max(a) < min_top_alpha) {
      next
    }

    ord <- order(a, decreasing = TRUE)
    k <- which(cumsum(a[ord]) >= coverage)[1]
    vars <- ord[seq_len(k)]
    p <- if (length(vars) == 1) {
      1
    } else {
      R_sub <- abs(R[vars, vars, drop = FALSE])
      min(R_sub[upper.tri(R_sub)])
    }

    if (is.finite(p) && p >= min_abs_corr) {
      nm <- paste0("L", ell)
      cs[[nm]] <- vars
      purity[nm] <- p
    }
  }

  list(cs = cs, purity = purity)
}

count_selected_cs <- function(alpha, R, coverage = 0.95,
                              min_abs_corr = 0.5,
                              min_top_alpha = 1e-4) {
  length(make_cs_from_alpha(
    alpha = alpha,
    R = R,
    coverage = coverage,
    min_abs_corr = min_abs_corr,
    min_top_alpha = min_top_alpha
  )$cs)
}

last_elbo <- function(fit) {
  tail(stats::na.omit(fit$elbo), 1)
}

sim_asn_eur <- function(chrom, start, end, n_gwas, n_ref, J, h2, seed) {
  if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
    python <- Sys.which("python3")
    if (nzchar(python)) {
      Sys.setenv(RETICULATE_PYTHON = python)
      reticulate::use_python(python, required = FALSE)
    }
  }

  stdpopsim <- reticulate::import("stdpopsim")
  msprime   <- reticulate::import("msprime")

  set.seed(seed)

  species <- stdpopsim$get_species("HomSap")
  contig  <- species$get_contig(chrom, left = as.integer(start), right = as.integer(end))

  demo_model <- species$get_demographic_model("OutOfAfrica_3G09")
  demography <- demo_model$model

  ts <- msprime$sim_ancestry(
    samples = reticulate::dict(CHB = as.integer(n_gwas), CEU = as.integer(n_ref)),
    demography = demography,
    recombination_rate = contig$recombination_map,
    random_seed = as.integer(seed)
  )

  ts <- msprime$sim_mutations(ts, rate = 1.29e-8, random_seed = as.integer(seed))

  chb_id <- NULL
  ceu_id <- NULL
  for (p in reticulate::iterate(ts$populations())) {
    if (!is.null(p$metadata) && p$metadata$name == "CHB") chb_id <- p$id
    if (!is.null(p$metadata) && p$metadata$name == "CEU") ceu_id <- p$id
  }

  samples_gwas <- ts$samples(population = as.integer(chb_id))
  samples_ref  <- ts$samples(population = as.integer(ceu_id))

  G_all <- ts$genotype_matrix()

  to_diploid <- function(hap_idx) {
    idx <- as.integer(hap_idx) + 1
    g_hap <- G_all[, idx, drop = FALSE]
    g_dip <- g_hap[, seq(1, ncol(g_hap), by = 2)] +
      g_hap[, seq(2, ncol(g_hap), by = 2)]
    t(g_dip)
  }

  G_gwas <- to_diploid(samples_gwas)
  G_ref  <- to_diploid(samples_ref)

  maf_gwas <- colMeans(G_gwas) / 2
  maf_ref  <- colMeans(G_ref) / 2

  common_mask <- (maf_gwas > 0.05) & (maf_gwas < 0.95) &
    (maf_ref > 0.05) & (maf_ref < 0.95)
  n_common <- sum(common_mask)
  if (n_common < J) {
    stop(paste0("Only ", n_common, " common variants found but J = ", J,
                ". Try a larger region."))
  }
  G_gwas <- G_gwas[, common_mask, drop = FALSE]
  G_ref  <- G_ref[, common_mask, drop = FALSE]

  J_all <- min(ncol(G_gwas), J)
  max_start <- J_all - J
  start_idx <- if (max_start > 0) sample(0:max_start, 1) else 0
  cols_to_keep <- (start_idx + 1):(start_idx + J)
  G  <- G_gwas[, cols_to_keep, drop = FALSE]
  G0 <- G_ref[, cols_to_keep, drop = FALSE]

  X <- sweep(G, 2, colMeans(G), "-")
  N <- nrow(X)

  G0_centered <- sweep(G0, 2, colMeans(G0), "-")

  R0 <- stats::cov(G0_centered)
  R0[is.na(R0)] <- 0

  R <- stats::cov(X)
  R[is.na(R)] <- 0

  num_causal <- 2
  first_causal_SNP <- sample(1:J, 1)
  second_causal_SNP <- which.min(R[first_causal_SNP, ]^2)
  causal_SNPs_loc <- c(first_causal_SNP, second_causal_SNP)

  beta_true <- rep(0, J)
  beta_true[causal_SNPs_loc] <- stats::runif(num_causal) + 0.01

  var_expl <- sum((X %*% beta_true)^2 / N)
  scale <- h2 / var_expl
  beta_true <- beta_true * sqrt(scale)

  y <- (X %*% beta_true) + sqrt(1 - h2) * stats::rnorm(N)
  y_std <- sqrt(sum((y - mean(y))^2) / N)
  y <- (y - mean(y)) / y_std

  v <- (t(X) %*% y) / N

  list(
    v = as.vector(v),
    R = R,
    R0 = R0,
    beta_true = beta_true,
    causal_SNPs_loc = causal_SNPs_loc
  )
}

standardize_diff_ancestry_sim <- function(res) {
  d <- 1 / sqrt(diag(res$R))
  res$v <- d * res$v
  res$R <- cov2cor(res$R)
  res$R0 <- cov2cor(res$R0)
  res
}

fit_relax_lambda_grid <- function(x, R0, N, L, beta_true, lambda_grid,
                                  sigma2_init = 0.3^2,
                                  max_iter = 3000) {
  fit_time <- system.time({
    fits <- lapply(lambda_grid, function(lambda) {
      fit_susie_relax(
        x = x,
        Rbar = make_rbar(R0, lambda = lambda),
        N = N,
        L = L,
        sigma2 = rep(sigma2_init, L),
        nu0_init = 2000,
        estimate_nu0 = TRUE,
        estimate_sigma2 = TRUE,
        warmup_iter = 5,
        max_iter = max_iter,
        tol = 1e-7,
        verbose = FALSE,
        nu0_bounds = c(10, 10000),
        sigma2_bounds = c(1e-12, 1),
        elbo_update_interval = 5,
        nu0_update_interval = 1,
        sigma2_update_interval = 1,
        nu0_tol = 1e-4
      )
    })
  })
  names(fits) <- sprintf("lambda=%g", lambda_grid)

  profile <- do.call(rbind, lapply(seq_along(lambda_grid), function(i) {
    fit <- fits[[i]]
    Rbar_i <- make_rbar(R0, lambda = lambda_grid[i])
    data.frame(
      method = "relax",
      lambda = lambda_grid[i],
      final_elbo = last_elbo(fit),
      n_iter = length(fit$elbo),
      nu0 = fit$nu0,
      selected_L = count_selected_cs(fit$alpha, R = Rbar_i),
      causal_pip = paste(sprintf("%.3f", fit$pip[beta_true != 0]), collapse = ", "),
      top = paste(apply(fit$alpha, 1, which.max), collapse = ", "),
      top_alpha = paste(sprintf("%.3f", apply(fit$alpha, 1, max)), collapse = ", ")
    )
  }))

  best_idx <- which.max(profile$final_elbo)
  list(
    best_fit = fits[[best_idx]],
    best_idx = best_idx,
    best_lambda = profile$lambda[best_idx],
    profile = profile,
    fit_seconds = unname(fit_time["elapsed"])
  )
}

fit_collapse_lambda_grid <- function(x, R0, N, L, beta_true, lambda_grid,
                                     nu0_grid,
                                     sigma2_init = 0.3^2,
                                     max_iter = 200) {
  fit_time <- system.time({
    fits <- lapply(lambda_grid, function(lambda) {
      fit_susie_collapse_r_shared(
        x = x,
        Rbar = make_rbar(R0, lambda = lambda),
        N = N,
        L = L,
        sigma2 = rep(sigma2_init, L),
        nu0_grid = nu0_grid,
        estimate_sigma2 = TRUE,
        max_iter = max_iter,
        tol = 1e-6,
        verbose = FALSE
      )
    })
  })
  names(fits) <- sprintf("lambda=%g", lambda_grid)

  profile <- do.call(rbind, lapply(seq_along(lambda_grid), function(i) {
    fit <- fits[[i]]
    Rbar_i <- make_rbar(R0, lambda = lambda_grid[i])
    data.frame(
      method = "collapse_shared",
      lambda = lambda_grid[i],
      final_elbo = fit$lower_bound,
      n_iter = length(fit$elbo),
      nu0 = fit$nu0,
      selected_L = count_selected_cs(fit$alpha, R = Rbar_i),
      causal_pip = paste(sprintf("%.3f", fit$pip[beta_true != 0]), collapse = ", "),
      top = paste(apply(fit$alpha, 1, which.max), collapse = ", "),
      top_alpha = paste(sprintf("%.3f", apply(fit$alpha, 1, max)), collapse = ", ")
    )
  }))

  best_idx <- which.max(profile$final_elbo)
  list(
    best_fit = fits[[best_idx]],
    best_idx = best_idx,
    best_lambda = profile$lambda[best_idx],
    profile = profile,
    fit_seconds = unname(fit_time["elapsed"])
  )
}

compare_one_seed <- function(seed, chrom = "chr22", start = 20e6, end = 21e6,
                             N = 25000L, N0 = 500L, J = 500L, h2 = 0.01,
                             L = 5,
                             lambda_grid = c(1e-4, 1e-3, 1e-2, 0.05, 0.1,
                                             0.15, 0.2, 0.3, 0.4, 0.5),
                             collapse_nu0_grid = c(10, 25, 50, 100, 200,
                                                   500, 1000, 2000, 5000),
                             relax_max_iter = 3000,
                             collapse_max_iter = 200,
                             verbose = TRUE) {
  if (verbose) {
    message(sprintf("seed %d: simulate", seed))
  }
  sim_time <- system.time({
    res <- sim_asn_eur(
      chrom = chrom,
      start = start,
      end = end,
      n_gwas = N,
      n_ref = N0,
      J = J,
      h2 = h2,
      seed = seed
    )
    res <- standardize_diff_ancestry_sim(res)
  })

  x <- res$v
  R0 <- res$R0
  beta_true <- res$beta_true
  causal <- which(beta_true != 0)

  if (verbose) {
    message(sprintf("seed %d: fit SuSiE-relax lambda grid", seed))
  }
  relax <- fit_relax_lambda_grid(
    x = x,
    R0 = R0,
    N = N,
    L = L,
    beta_true = beta_true,
    lambda_grid = lambda_grid,
    max_iter = relax_max_iter
  )

  if (verbose) {
    message(sprintf("seed %d: fit collapsed-R lambda grid", seed))
  }
  collapse <- fit_collapse_lambda_grid(
    x = x,
    R0 = R0,
    N = N,
    L = L,
    beta_true = beta_true,
    lambda_grid = lambda_grid,
    nu0_grid = collapse_nu0_grid,
    max_iter = collapse_max_iter
  )

  best_rows <- rbind(
    transform(relax$profile[relax$best_idx, ], seed = seed),
    transform(collapse$profile[collapse$best_idx, ], seed = seed)
  )
  best_rows$causal <- paste(causal, collapse = ", ")
  best_rows$sim_seconds <- unname(sim_time["elapsed"])
  best_rows$fit_seconds <- c(relax$fit_seconds, collapse$fit_seconds)

  list(
    seed = seed,
    best = best_rows,
    profiles = rbind(
      transform(relax$profile, seed = seed),
      transform(collapse$profile, seed = seed)
    )
  )
}

run_diff_ancestry_comparison <- function(
    seeds = 1:25,
    output_prefix = "analysis/compare-collapse-vs-relax-diff-ancestry",
    ...) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)

  runs <- vector("list", length(seeds))
  names(runs) <- paste0("seed", seeds)

  for (i in seq_along(seeds)) {
    runs[[i]] <- compare_one_seed(seed = seeds[i], ...)
    partial_best <- do.call(rbind, lapply(runs[seq_len(i)], `[[`, "best"))
    partial_profiles <- do.call(rbind, lapply(runs[seq_len(i)], `[[`, "profiles"))
    utils::write.csv(partial_best, paste0(output_prefix, "-best.csv"), row.names = FALSE)
    utils::write.csv(partial_profiles, paste0(output_prefix, "-profiles.csv"), row.names = FALSE)
    saveRDS(runs, paste0(output_prefix, "-runs.rds"))
  }

  best <- do.call(rbind, lapply(runs, `[[`, "best"))
  profiles <- do.call(rbind, lapply(runs, `[[`, "profiles"))
  selected_L_distribution <- as.data.frame.matrix(table(best$method, best$selected_L))
  lambda_distribution <- as.data.frame.matrix(table(best$method, best$lambda))
  nu0_distribution <- as.data.frame.matrix(table(best$method, best$nu0))

  out <- list(
    best = best,
    profiles = profiles,
    selected_L_distribution = selected_L_distribution,
    lambda_distribution = lambda_distribution,
    nu0_distribution = nu0_distribution,
    runs = runs
  )
  saveRDS(out, paste0(output_prefix, "-summary.rds"))
  out
}

print_comparison_summary <- function(cmp) {
  cat("\nSelected L distribution:\n")
  print(cmp$selected_L_distribution)
  cat("\nBest lambda distribution:\n")
  print(cmp$lambda_distribution)
  cat("\nBest nu0 distribution:\n")
  print(cmp$nu0_distribution)
  cat("\nPer-seed best fits:\n")
  print(cmp$best)
  invisible(cmp)
}

if (identical(Sys.getenv("RUN_DIFF_ANCESTRY_COMPARISON"), "true")) {
  cmp <- run_diff_ancestry_comparison(seeds = 1:25)
  print_comparison_summary(cmp)
}
