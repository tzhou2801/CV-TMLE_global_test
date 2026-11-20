# =============================================================================
# Simulation Study 2 (Section 6.2 and Appendix A) for:
#   "A CV-TMLE global test approach to improve power in rare disease clinical
#    studies with multiple-component endpoints" (arXiv:2605.02851)
# =============================================================================
#
# This script implements a simulation study comparing FOUR global null hypothesis
# testing methods for TWO correlated outcomes, generated to mimic a realistic
# rare-disease trial described in Section 14.1 / Table 5 of the ELAPRASE
# (idursulfase) prescribing information:
#   - Two arms (placebo vs active), fully observed
#   - Baseline covariates: age, region, baseline disease severity
#   - Baseline outcomes: 6MWT and %-predicted FVC with entry criterion FVC < 80%
#   - Week 53 outcomes generated via (baseline + change)
#   - Analysis endpoints include both:
#       * Continuous changes: d6MWT = change in 6MWT, dFVC = change in FVC
#       * Rank changes:       Y1 = rank(d6MWT), Y2 = rank(dFVC)
#
# METHODS COMPARED:
# 1. Holm correction applied to TMLE p-values for (d6MWT, dFVC)
# 2. Hochberg procedure applied to TMLE p-values for (d6MWT, dFVC)
# 3. O'Brien global rank-sum test on rank outcomes (Y1, Y2) (UNADJUSTED; original O'Brien)
# 4. Stabilized CV-TMLE on (d6MWT, dFVC) with data-adaptive weights (includes shrinkage step)
#
# STABILIZATION P-VALUE:
# The stabilization p-value in the CV-TMLE method uses a Monte Carlo
# approximation to the null distribution of the supremum statistic
# T*_v = sup_{alpha} T_alpha (Section 4.5).
#
# Author: Tianyue Zhou, tianyue_zhou@berkeley.edu
# Date: March 2026
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(MASS)          # mvrnorm
  library(tmle)
  library(SuperLearner)
  library(foreach)
  library(doParallel)
})

# =============================================================================
# 1. REALISTIC DATA GENERATING PROCESS (DGP)
# =============================================================================

# Generate one 2-arm trial dataset with:
#   - Covariates: age (5-31), region (3-level factor), severity (3-level factor)
#   - Baseline outcomes (B_6MWT, B_FVC) with truncation FVC < 80
#   - Week53 outcomes via change distribution (Delta_6MWT, Delta_FVC)
#   - Analysis outcomes: continuous changes (d6MWT, dFVC) and ranks (Y1, Y2)
#
# Scenario:
#   - "alt": Table-5-like treatment effect (active better mean changes)
#   - "weak_alt": weak joint alternative (both endpoints have modest effects; no single endpoint dominates)
#   - "null": global null (active mean changes = placebo mean changes)
#
simulate_rank_change_trial <- function(
  n_total = 60,
  scenario = c("alt", "weak_alt", "null"),
  seed = NULL,

  # ---- Targets from Table 5 (marginal, by arm) ----
  mu_B = c(392.5, 55.45),        # baseline mean (6MWT, FVC)
  sd_B = c(107, 14),
  rho_B = 0.30,

  mu_D_placebo = c(7, 0.8),
  mu_D_active_alt = c(44, 3.4),
  mu_D_active_weak = NULL,          # weak joint alt mean changes; if NULL defaults to (mu_D_placebo + c(25, 2.6))
  sd_D_placebo = c(54, 9.6),
  sd_D_active  = c(70, 10.0),
  rho_D = 0.25,                  # assumed (not reported in the label)

  # Under the global-null scenario, you may want *identical* change distributions
  # in both arms (not just equal means). Setting null_equal_variance=TRUE makes the
  # active-arm change SDs equal the placebo SDs when scenario=="null".
  # This improves type-I error calibration for distribution-based tests (e.g., O'Brien).
  null_equal_variance = TRUE,

  # ---- Covariate generation (from label: ages 5-31) ----
  age_mean = 15,
  age_sd = 6,
  age_bounds = c(5, 31),

  region_levels = c("NA", "EU", "Other"),
  region_probs  = c(0.4, 0.4, 0.2),

  # ---- Disease severity category (3-level observed baseline covariate) ----
  # We generate sev_cat directly (no latent continuous score). By default, we choose
  # probabilities implied by cutting a N(0,1) at sev_cut, to preserve the old mix.
  sev_cut = c(-0.5, 0.5),        # only used to define default sev_probs when sev_probs=NULL
  sev_probs = NULL,              # optional override: c(p_mild, p_moderate, p_severe)

  # ---- Covariate effects on BASELINE outcomes (stronger; depends on sev_cat) ----
  # X = (age10, sev_mod_c, sev_sev_c, EU_c, Other_c)
  #   age10      = (age - mean(age))/10
  #   sev_mod_c  = I(sev_cat=="moderate") - mean()
  #   sev_sev_c  = I(sev_cat=="severe")   - mean()
  #   EU_c/Other_c are centered region indicators.
  #
  # Coefficients are interpretable as differences vs "mild" because indicators are centered.
  # Example: baseline 6MWT moderate is ~60m lower than mild, severe ~120m lower than mild.
  M_B = matrix(c(
    -15,  -60, -120,  -10,   5,     # baseline 6MWT shifts
    -1.5,  -6,  -12, -1.5, 0.75     # baseline FVC shifts
  ), nrow = 2, byrow = TRUE),

  # ---- Covariate effects on CHANGE outcomes (depends on sev_cat) ----
  # Example: change 6MWT moderate is ~18m lower than mild, severe ~40m lower than mild.
  M_D = matrix(c(
     -8,  -18,  -40,   -4,   2,     # change in 6MWT
     -1,   -5,  -10, -0.8, 0.5      # change in FVC
  ), nrow = 2, byrow = TRUE),

  # Optional global multipliers (easy knob for sensitivity analysis)
  M_B_mult = 1,
  M_D_mult = 1,

  # ---- Plausibility / inclusion bounds ----
  bounds_6mwt = c(50, 650),
  bounds_fvc  = c(20, 80)         # enforce 20 <= FVC < 80
) {
  scenario <- match.arg(scenario)
  if (!is.null(seed)) set.seed(seed)

  # Default weak alternative (anchored to placebo mean changes).
  # If you change mu_D_placebo, this keeps weak-alt differences anchored to placebo.
  if (is.null(mu_D_active_weak)) {
    mu_D_active_weak <- mu_D_placebo + c(25, 2.6)   # +25m 6MWT, +2.6pp FVC vs placebo
  }

  # Balanced randomization (as close as possible)
  A <- sample(rep(c(0, 1), length.out = n_total))  # 0=placebo, 1=active

  # ---- helper: truncated normal for age
  rtruncnorm1 <- function(mu, sd, lo, hi) {
    repeat {
      x <- rnorm(1, mean = mu, sd = sd)
      if (x >= lo && x <= hi) return(x)
    }
  }

  # ---- covariates
  age <- vapply(
    1:n_total,
    function(i) rtruncnorm1(age_mean, age_sd, age_bounds[1], age_bounds[2]),
    numeric(1)
  )

  region <- factor(
    sample(region_levels, n_total, replace = TRUE, prob = region_probs),
    levels = region_levels
  )

  # Severity category: generate directly (no latent score)
  if (is.null(sev_probs)) {
    sev_probs <- c(
      pnorm(sev_cut[1]),
      pnorm(sev_cut[2]) - pnorm(sev_cut[1]),
      1 - pnorm(sev_cut[2])
    )
  }
  sev_probs <- sev_probs / sum(sev_probs)

  sev_cat <- factor(
    sample(c("mild", "moderate", "severe"), n_total, replace = TRUE, prob = sev_probs),
    levels = c("mild", "moderate", "severe")
  )

  # ---- centered design columns X = (age10, sev_mod_c, sev_sev_c, EU_c, Other_c)
  # IMPORTANT: we center using *fixed (population)* means/probabilities rather than
  # sample means, to avoid a subtle across-subject dependence induced by centering on
  # the realized sample mean within each dataset.
  truncnorm_mean <- function(mu, sd, lo, hi) {
    a <- (lo - mu) / sd
    b <- (hi - mu) / sd
    Z <- pnorm(b) - pnorm(a)
    mu + sd * (dnorm(a) - dnorm(b)) / Z
  }

  age_center <- truncnorm_mean(age_mean, age_sd, age_bounds[1], age_bounds[2])
  age10 <- (age - age_center) / 10

  sev_mod <- as.numeric(sev_cat == "moderate")
  sev_sev <- as.numeric(sev_cat == "severe")
  sev_mod_c <- sev_mod - sev_probs[2]
  sev_sev_c <- sev_sev - sev_probs[3]

  EU    <- as.numeric(region == "EU")
  Other <- as.numeric(region == "Other")
  EU_c    <- EU    - region_probs[2]
  Other_c <- Other - region_probs[3]

  X <- cbind(
    age10 = age10,
    sev_mod_c = sev_mod_c,
    sev_sev_c = sev_sev_c,
    EU_c = EU_c,
    Other_c = Other_c
  )

  # Apply optional multipliers (sensitivity knob)
  M_B <- M_B_mult * M_B
  M_D <- M_D_mult * M_D

  # ---- target covariance matrices helper
  Sigma_target <- function(sd_vec, rho) {
    matrix(c(sd_vec[1]^2, rho * sd_vec[1] * sd_vec[2],
             rho * sd_vec[1] * sd_vec[2], sd_vec[2]^2),
           2, 2)
  }
  Sigma_B_target  <- Sigma_target(sd_B, rho_B)
  Sigma_D0_target <- Sigma_target(sd_D_placebo, rho_D)

  sd_D_active_use <- sd_D_active
  if (scenario == "null" && isTRUE(null_equal_variance)) {
    sd_D_active_use <- sd_D_placebo
  }
  Sigma_D1_target <- Sigma_target(sd_D_active_use, rho_D)

  # ---- compute residual covariance to keep (approximately) desired marginals
  # Sigma_resid = Sigma_target - M * Var(X) * M^T  (shrink M if needed)
  make_resid <- function(Sigma_target_mat, M, Xmat, min_eig = 1e-8) {
    # NOTE: cov() always centers internally; Xmat is population-centered (expected mean 0).
    Sx <- stats::cov(Xmat)
    M_use <- M
    repeat {
      Sigma_cov <- M_use %*% Sx %*% t(M_use)
      Sigma_res <- Sigma_target_mat - Sigma_cov
      Sigma_res <- (Sigma_res + t(Sigma_res)) / 2
      ev <- eigen(Sigma_res, symmetric = TRUE)$values
      if (min(ev) > min_eig) return(list(Sigma_resid = Sigma_res, M_use = M_use))
      M_use <- 0.9 * M_use  # shrink covariate effects until PSD
    }
  }

  base_res <- make_resid(Sigma_B_target,  M_B, X)
  Sigma_B_resid <- base_res$Sigma_resid
  M_B_use <- base_res$M_use

  chg0_res <- make_resid(Sigma_D0_target, M_D, X)
  Sigma_D0_resid <- chg0_res$Sigma_resid
  M_D_use <- chg0_res$M_use

  chg1_res <- make_resid(Sigma_D1_target, M_D_use, X)
  Sigma_D1_resid <- chg1_res$Sigma_resid

  # ---- helper: truncated BVN for baseline given subject-specific mean
  rmvnorm_trunc1 <- function(mu_i, Sigma, bounds_6mwt, bounds_fvc, max_iter = 2e5) {
    for (iter in 1:max_iter) {
      x <- MASS::mvrnorm(1, mu = mu_i, Sigma = Sigma)
      ok <- (x[1] >= bounds_6mwt[1] && x[1] <= bounds_6mwt[2] &&
             x[2] >= bounds_fvc[1]  && x[2] <  bounds_fvc[2])
      if (ok) return(x)
    }
    stop("Baseline truncation failed; loosen bounds or reduce covariate effects (M_B).")
  }

  # ---- baseline outcomes
  B <- matrix(NA_real_, nrow = n_total, ncol = 2)
  for (i in 1:n_total) {
    mu_i <- as.numeric(mu_B + M_B_use %*% X[i, ])
    B[i, ] <- rmvnorm_trunc1(mu_i, Sigma_B_resid, bounds_6mwt, bounds_fvc)
  }
  colnames(B) <- c("B_6MWT", "B_FVC")

  # ---- scenario-specific mean change for active arm
  mu_D_active <- switch(
    scenario,
    null     = mu_D_placebo,
    alt      = mu_D_active_alt,
    weak_alt = mu_D_active_weak
  )

  # ---- changes (fully observed)
  Delta <- matrix(NA_real_, nrow = n_total, ncol = 2)
  for (i in 1:n_total) {
    if (A[i] == 0) {
      muD_i <- as.numeric(mu_D_placebo + M_D_use %*% X[i, ])
      Delta[i, ] <- MASS::mvrnorm(1, mu = muD_i, Sigma = Sigma_D0_resid)
    } else {
      muD_i <- as.numeric(mu_D_active  + M_D_use %*% X[i, ])
      Delta[i, ] <- MASS::mvrnorm(1, mu = muD_i, Sigma = Sigma_D1_resid)
    }
  }
  colnames(Delta) <- c("D_6MWT", "D_FVC")

  # Week 53 outcomes (not strictly needed for testing, but kept for realism)
  Y53 <- B + Delta
  colnames(Y53) <- c("Y53_6MWT", "Y53_FVC")

  # ---- analysis outcomes: ranks of changes (higher = better improvement)
  Y1_rank <- rank(Delta[, 1], ties.method = "average")
  Y2_rank <- rank(Delta[, 2], ties.method = "average")

  data.frame(
    A = A,
    age = age,
    region = region,
    sev_cat = sev_cat,
    B_6MWT = B[, 1],
    B_FVC = B[, 2],
    Y53_6MWT = Y53[, 1],
    Y53_FVC = Y53[, 2],
    d6MWT = Delta[, 1],
    dFVC = Delta[, 2],
    Y1 = Y1_rank,   # rank(change 6MWT)
    Y2 = Y2_rank    # rank(change FVC)
  )
}



# =============================================================================
# 2. STATISTICAL TEST FUNCTIONS
# =============================================================================

# -------------------------
# Methods 1 & 2: TMLE with Holm and Hochberg corrections
# -------------------------
tmle_multiplicity_test <- function(data, treatment_var, outcome1_var, outcome2_var,
                                   covariate_vars = c("age", "region", "sev_cat"),
                                   alpha = 0.025,
                                   one_sided = TRUE,
                                   Q.SL.library = c("SL.glm"),
                                   verbose = FALSE,
                                   seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(data)
  df_t_use <- max(1, n - 4)

  Y1 <- data[[outcome1_var]]
  Y2 <- data[[outcome2_var]]
  A <- data[[treatment_var]]
  W <- data[, covariate_vars, drop = FALSE]

  # Known propensity score vector (g0 = 0.5 for RCT)
  g0 <- rep(0.5, n)

  tmle1 <- tmle(Y = Y1, A = A, W = W,
                Q.SL.library = Q.SL.library,
                g1W = g0,
                family = "gaussian",
                verbose = verbose)

  tmle2 <- tmle(Y = Y2, A = A, W = W,
                Q.SL.library = Q.SL.library,
                g1W = g0,
                family = "gaussian",
                verbose = verbose)

  se1 <- sqrt(tmle1$estimates$ATE$var.psi)
  se2 <- sqrt(tmle2$estimates$ATE$var.psi)
  z1 <- tmle1$estimates$ATE$psi / se1
  z2 <- tmle2$estimates$ATE$psi / se2

  if (one_sided) {
    p1 <- pt(z1, df = df_t_use, lower.tail = FALSE)
    p2 <- pt(z2, df = df_t_use, lower.tail = FALSE)
  } else {
    p1 <- 2 * pt(abs(z1), df = df_t_use, lower.tail = FALSE)
    p2 <- 2 * pt(abs(z2), df = df_t_use, lower.tail = FALSE)
  }

  p_values <- c(p1, p2)

  # Holm (m=2: equivalent to Bonferroni)
  holm_reject <- (p1 < alpha / 2) || (p2 < alpha / 2)

  # Hochberg step-up
  m <- length(p_values)
  p_order <- order(p_values, decreasing = TRUE)
  p_sorted <- p_values[p_order]
  hochberg_reject <- FALSE
  for (i in 1:m) {
    if (p_sorted[i] < alpha / i) {
      hochberg_reject <- TRUE
      break
    }
  }

  list(
    holm_reject = holm_reject,
    hochberg_reject = hochberg_reject,
    p_values = p_values,
    estimates = c(tmle1$estimates$ATE$psi, tmle2$estimates$ATE$psi),
    standard_errors = c(se1, se2),
    z_values = c(z1, z2),
    df_t = df_t_use
  )
}

# -------------------------
# Method 3: O'Brien global rank-sum test (permutation inference)
# -------------------------
# Inputs are the *rank outcomes* (Y1,Y2). O'Brien score is the sum of endpoint ranks.
obrien_rank_sum_test <- function(data, treatment_var, outcome_vars,
                                 alpha = 0.025,
                                 one_sided = TRUE,
                                 n_perm = 5000,
                                 seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  A <- data[[treatment_var]]
  stopifnot(all(outcome_vars %in% names(data)))

  # O'Brien score: sum of endpoint ranks (higher = better)
  score <- rowSums(data[, outcome_vars, drop = FALSE])
  n1 <- sum(A == 1)
  n0 <- sum(A == 0)

  # Statistic: difference in mean scores (active - placebo)
  stat_obs <- mean(score[A == 1]) - mean(score[A == 0])

  # Permutation (randomization) inference:
  #   permute treatment labels while holding n1 fixed.
  #   Use +1 correction for an exact (finite-sample) p-value under randomization.
  n <- length(score)
  total_score <- sum(score)

  # Fast permutation: only need sum of treated scores for each permutation
  perm_sum_treated <- replicate(n_perm, sum(score[sample.int(n, n1, replace = FALSE)]))
  stat_perm <- perm_sum_treated / n1 - (total_score - perm_sum_treated) / n0

  if (one_sided) {
    p_value <- (1 + sum(stat_perm >= stat_obs)) / (n_perm + 1)
  } else {
    p_value <- (1 + sum(abs(stat_perm) >= abs(stat_obs))) / (n_perm + 1)
  }

  list(
    reject_h0 = (p_value < alpha),
    p_value = p_value,
    stat_obs = stat_obs,
    n_perm = n_perm
  )
}

# -------------------------
# Monte Carlo supremum p-value helper
# -------------------------

#' Compute Monte Carlo p-value for the supremum statistic T* = sup_{alpha} T_alpha
#' under the null distribution Z ~ N(0, Sigma).
#'
#' @param Sigma_mat 2x2 covariance matrix of influence curves (from training data)
#' @param observed_T The observed supremum statistic T*_v
#' @param B Number of Monte Carlo draws (default 5000)
#' @return Monte Carlo p-value P(T*_mc > observed_T)
mc_supremum_pvalue <- function(Sigma_mat, observed_T, B = 5000) {
  Z_draws <- mvrnorm(n = B, mu = c(0, 0), Sigma = Sigma_mat)

  s11 <- Sigma_mat[1, 1]
  s12 <- Sigma_mat[1, 2]
  s22 <- Sigma_mat[2, 2]

  T_at_0 <- Z_draws[, 2] / sqrt(s22)
  T_at_1 <- Z_draws[, 1] / sqrt(s11)
  T_star_mc <- pmax(T_at_0, T_at_1)

  for (b in 1:B) {
    z1 <- Z_draws[b, 1]
    z2 <- Z_draws[b, 2]
    opt <- optimize(
      f = function(a) {
        num <- a * z1 + (1 - a) * z2
        denom <- sqrt(a^2 * s11 + 2 * a * (1 - a) * s12 + (1 - a)^2 * s22)
        return(-num / denom)
      },
      interval = c(0, 1)
    )
    T_interior <- -opt$objective
    if (T_interior > T_star_mc[b]) T_star_mc[b] <- T_interior
  }

  p_val <- mean(T_star_mc > observed_T)
  return(p_val)
}

# -------------------------
# Method 4: Stabilized CV-TMLE with data-adaptive weights
# -------------------------
cv_tmle_test <- function(data, treatment_var, outcome1_var, outcome2_var,
                         covariate_vars = c("age", "region", "sev_cat"),
                         alpha = 0.025,
                         V = 10,
                         Q.SL.library = c("SL.glm"),
                         verbose = FALSE,
                         seed = NULL,

                         # Shrinkage (stabilization) parameters (Sections 4.5-4.6)
                         alpha_ref = c(0.5, 0.5),
                         C = 2,                         # lambda = clamp(C*log(n)*p_n, 0, 1)
                         df_t = 36                     # final critical value uses t_df (Logan-Tamhane)
) {

  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(data)
  Y1 <- data[[outcome1_var]]
  Y2 <- data[[outcome2_var]]
  A <- data[[treatment_var]]
  W <- data[, covariate_vars, drop = FALSE]

  # Degrees of freedom for t-distribution-based calibration
  # Avoid df > (n-2) if n is very small (not expected here, but makes function robust).
  df_t_use <- min(df_t, n - 2)

  # Create folds
  folds <- sample(rep(1:V, length.out = n))

  # Storage
  psi_star_v <- numeric(V)
  IC_all <- numeric(n)

  alpha_v_matrix <- matrix(NA, nrow = V, ncol = 2)     # stabilized weights per fold
  alpha_adapt_matrix <- matrix(NA, nrow = V, ncol = 2) # adaptive weights (pre-shrink)
  p_shrink_v <- rep(NA_real_, V)                       # p_n,v
  lambda_v <- rep(NA_real_, V)                         # shrinkage coefficient

  # Helper: clamp
  clamp01 <- function(x) max(0, min(1, x))

  # Process each fold
  for (v in 1:V) {
    if (verbose) cat("Processing fold", v, "of", V, "\n")

    train_idx <- which(folds != v)
    valid_idx <- which(folds == v)

    train_data <- data[train_idx, ]
    valid_data <- data[valid_idx, ]

    # Known propensity scores for an RCT
    g0_train <- rep(0.5, nrow(train_data))

    # Step 1: individual TMLEs on training data (continuous changes; gaussian working model)
    tmle1_train <- tmle(
      Y = train_data[[outcome1_var]],
      A = train_data[[treatment_var]],
      W = train_data[, covariate_vars, drop = FALSE],
      Q.SL.library = Q.SL.library,
      g1W = g0_train,
      family = "gaussian",
      verbose = FALSE
    )

    tmle2_train <- tmle(
      Y = train_data[[outcome2_var]],
      A = train_data[[treatment_var]],
      W = train_data[, covariate_vars, drop = FALSE],
      Q.SL.library = Q.SL.library,
      g1W = g0_train,
      family = "gaussian",
      verbose = FALSE
    )

    psi1_star <- tmle1_train$estimates$ATE$psi
    psi2_star <- tmle2_train$estimates$ATE$psi

    IC1_train <- tmle1_train$estimates$IC$IC.ATE
    IC2_train <- tmle2_train$estimates$IC$IC.ATE

    sigma11 <- var(IC1_train)
    sigma22 <- var(IC2_train)
    sigma12 <- cov(IC1_train, IC2_train)

    # Step 1a: L-BFGS-B optimization for optimal alpha (maximize test statistic)
    Sigma_mat <- matrix(c(sigma11, sigma12, sigma12, sigma22), 2, 2)

    objective_fn <- function(alpha1) {
      alpha2 <- 1 - alpha1
      psi_combined <- alpha1 * psi1_star + alpha2 * psi2_star
      alpha_full <- c(alpha1, alpha2)
      var_combined <- as.numeric(t(alpha_full) %*% Sigma_mat %*% alpha_full)
      if (var_combined > 1e-10) {
        return(-psi_combined / sqrt(var_combined))
      } else {
        return(1e10)
      }
    }

    gradient_fn <- function(alpha1) {
      alpha2 <- 1 - alpha1
      alpha_full <- c(alpha1, alpha2)
      psi_combined <- alpha1 * psi1_star + alpha2 * psi2_star
      var_combined <- as.numeric(t(alpha_full) %*% Sigma_mat %*% alpha_full)
      if (var_combined < 1e-10) return(0)

      sigma_combined <- sqrt(var_combined)
      Sigma_alpha <- as.numeric(Sigma_mat %*% alpha_full)
      dt_dalpha <- c(psi1_star, psi2_star) / sigma_combined -
        psi_combined * Sigma_alpha / (sigma_combined^3)

      # Chain rule: d/d(alpha1) = dt/d(alpha1) - dt/d(alpha2)
      return(-(dt_dalpha[1] - dt_dalpha[2]))
    }

    optim_result <- tryCatch({
      optim(par = 0.5, fn = objective_fn, gr = gradient_fn,
            method = "L-BFGS-B", lower = 0, upper = 1,
            control = list(maxit = 200, factr = 1e7))
    }, error = function(e) {
      list(par = 0.5, value = -objective_fn(0.5))
    })

    alpha1_opt <- optim_result$par
    alpha_adapt <- c(max(0, alpha1_opt), max(0, 1 - alpha1_opt))
    alpha_adapt <- alpha_adapt / sum(alpha_adapt)
    T_star <- -optim_result$value

    # Step 1b: shrinkage / stabilization
    p_n <- mc_supremum_pvalue(Sigma_mat, sqrt(nrow(train_data)) * T_star, B = 5000)
    lam <- clamp01(C * log(n) * p_n)

    alpha_stab <- (1 - lam) * alpha_adapt + lam * alpha_ref

    alpha_adapt_matrix[v, ] <- alpha_adapt
    alpha_v_matrix[v, ] <- alpha_stab
    p_shrink_v[v] <- p_n
    lambda_v[v] <- lam

    # Step 2: apply stabilized weights to validation data
    Y_bar_train <- alpha_stab[1] * train_data[[outcome1_var]] +
      alpha_stab[2] * train_data[[outcome2_var]]
    Y_bar_valid <- alpha_stab[1] * valid_data[[outcome1_var]] +
      alpha_stab[2] * valid_data[[outcome2_var]]

    # Outcome model for Y_bar (simple linear model for speed / comparability)
    train_model_df <- data.frame(
      Y_bar = Y_bar_train,
      A = train_data[[treatment_var]],
      train_data[, covariate_vars, drop = FALSE]
    )

    lm_formula <- as.formula(paste("Y_bar ~ A +", paste(covariate_vars, collapse = " + ")))
    lm_composite <- lm(lm_formula, data = train_model_df)

    valid_model_df <- data.frame(
      A = valid_data[[treatment_var]],
      valid_data[, covariate_vars, drop = FALSE]
    )

    Q_offset_v <- predict(lm_composite, newdata = valid_model_df)

    valid_model_df_1 <- valid_model_df
    valid_model_df_1$A <- 1
    Q_1_v <- predict(lm_composite, newdata = valid_model_df_1)

    valid_model_df_0 <- valid_model_df
    valid_model_df_0$A <- 0
    Q_0_v <- predict(lm_composite, newdata = valid_model_df_0)

    # Stack for pooled targeting
    if (v == 1) {
      Y_bar_stacked <- Y_bar_valid
      Q_offset_stacked <- Q_offset_v
      Q_1_stacked <- Q_1_v
      Q_0_stacked <- Q_0_v
      A_stacked <- valid_data[[treatment_var]]
      fold_id <- rep(v, length(valid_idx))
    } else {
      Y_bar_stacked <- c(Y_bar_stacked, Y_bar_valid)
      Q_offset_stacked <- c(Q_offset_stacked, Q_offset_v)
      Q_1_stacked <- c(Q_1_stacked, Q_1_v)
      Q_0_stacked <- c(Q_0_stacked, Q_0_v)
      A_stacked <- c(A_stacked, valid_data[[treatment_var]])
      fold_id <- c(fold_id, rep(v, length(valid_idx)))
    }
  }

  # Pooled targeting (known g0=0.5)
  H_stacked <- (2 * A_stacked - 1) / 0.5   # H=2 if A=1, H=-2 if A=0

  pooled_target <- lm(Y_bar_stacked ~ -1 + H_stacked, offset = Q_offset_stacked)
  epsilon <- coef(pooled_target)

  # Updated predictions
  Q_star_stacked <- Q_offset_stacked + epsilon * H_stacked
  Q_star_1_stacked <- Q_1_stacked + epsilon * 2       # H=2 when A=1
  Q_star_0_stacked <- Q_0_stacked + epsilon * (-2)    # H=-2 when A=0

  # Fold-specific psi and IC
  for (v in 1:V) {
    fold_mask <- (fold_id == v)

    Q_star_1_fold <- Q_star_1_stacked[fold_mask]
    Q_star_0_fold <- Q_star_0_stacked[fold_mask]
    Q_star_fold <- Q_star_stacked[fold_mask]
    Y_bar_fold <- Y_bar_stacked[fold_mask]
    H_fold <- H_stacked[fold_mask]

    psi_star_v[v] <- mean(Q_star_1_fold - Q_star_0_fold)

    IC_fold <- H_fold * (Y_bar_fold - Q_star_fold) +
      (Q_star_1_fold - Q_star_0_fold) - psi_star_v[v]

    IC_all[which(folds == v)] <- IC_fold
  }

  # Final CV-TMLE estimate and inference
  psi_cv_tmle <- mean(psi_star_v)
  sigma_cv <- sqrt(var(IC_all))

  t_stat <- sqrt(n) * psi_cv_tmle / sigma_cv
  p_value <- pt(t_stat, df = df_t_use, lower.tail = FALSE)

  alpha_avg <- colMeans(alpha_v_matrix)
  names(alpha_avg) <- c("alpha1", "alpha2")

  list(
    reject_h0 = (t_stat > qt(1 - alpha, df = df_t_use)),
    p_value = p_value,
    psi_cv_tmle = psi_cv_tmle,
    sigma_cv = sigma_cv,
    t_stat = t_stat,
    df_t = df_t_use,
    alpha_weights = alpha_avg,
    alpha_weights_matrix = alpha_v_matrix,
    alpha_adaptive_matrix = alpha_adapt_matrix,
    p_shrink = p_shrink_v,
    lambda_shrink = lambda_v
  )
}

# =============================================================================
# 3. SIMULATION PIPELINE (three scenarios: distribution-null global null + Table5 alternative + weak alternative)
# =============================================================================

run_four_method_simulation <- function(
  n_total = 60,
  scenarios = list(
    list(name = "global_null", scenario = "null"),
    list(name = "table5_alt",  scenario = "alt"),
    list(name = "weak_alt",    scenario = "weak_alt")
  ),
  n_reps = 1000,
  alpha = 0.025,
  covariate_vars = c("age", "region", "sev_cat"),
  tmle_libraries = list(Q.SL.library = c("SL.glm")),
  cv_tmle_params = list(V = 10, C = 2, df_t = 36),
  seed = 202701,
  verbose = TRUE,
  tmle_verbose = FALSE,
  n_cores = parallel::detectCores() - 1
) {

  set.seed(seed)

  if (verbose) {
    cat("Starting FOUR-METHOD simulation study...\n")
    cat("Methods: TMLE-Holm, TMLE-Hochberg, O'Brien(Unadj), Stabilized CV-TMLE\n")
    cat("Scenarios:", paste(vapply(scenarios, `[[`, "", "name"), collapse = ", "), "\n")
    cat("Replications:", n_reps, " | n_total:", n_total, " | alpha:", alpha, "\n")
    cat("Using", n_cores, "cores for parallel execution\n\n")
  }

  # ---- worker: run ALL scenarios for a single replication index
  one_replication <- function(rep,
                              scenarios, n_total, alpha,
                              covariate_vars, tmle_libraries, cv_tmle_params,
                              seed, tmle_verbose) {

    out_list <- vector("list", length(scenarios))

    for (s_idx in seq_along(scenarios)) {
      sc <- scenarios[[s_idx]]

      # Generate one dataset under this scenario (contains continuous changes and rank versions)
      sim_data <- simulate_rank_change_trial(
        n_total = n_total,
        scenario = sc$scenario,
        seed = seed + rep * 1000 + s_idx
      )

      # Methods 1 & 2: TMLE with Holm and Hochberg corrections
      tmle_result <- tmle_multiplicity_test(
        sim_data, "A", "d6MWT", "dFVC",
        covariate_vars = covariate_vars,
        alpha = alpha,
        one_sided = TRUE,
        Q.SL.library = tmle_libraries$Q.SL.library,
        verbose = tmle_verbose,
        seed = seed + rep * 1000 + 100 + s_idx
      )

      # Method 3
      obrien_result <- obrien_rank_sum_test(
        data = sim_data,
        treatment_var = "A",
        outcome_vars = c("Y1", "Y2"),
        alpha = alpha,
        one_sided = TRUE,
        seed = seed + rep * 1000 + 500 + s_idx
      )

      # Method 4 (stabilized CV-TMLE)
      cv_tmle_result <- do.call(
        cv_tmle_test,
        c(list(
          data = sim_data,
          treatment_var = "A",
          outcome1_var = "d6MWT",
          outcome2_var = "dFVC",
          covariate_vars = covariate_vars,
          alpha = alpha,
          Q.SL.library = tmle_libraries$Q.SL.library,
          verbose = FALSE,
          seed = seed + rep * 1000 + 999 + s_idx
        ), cv_tmle_params)
      )

      out_list[[s_idx]] <- data.frame(
        replication = rep,
        scenario_name = sc$name,
        n_total = n_total,
        observed_correlation = cor(sim_data$d6MWT, sim_data$dFVC),
        observed_correlation_rank = cor(sim_data$Y1, sim_data$Y2),

        # Rejections
        tmle_holm_reject = tmle_result$holm_reject,
        tmle_hochberg_reject = tmle_result$hochberg_reject,
        obrien_reject = obrien_result$reject_h0,
        cv_tmle_reject = cv_tmle_result$reject_h0,

        # P-values (TMLE methods use continuous changes)
        p_d6MWT = tmle_result$p_values[1],
        p_dFVC  = tmle_result$p_values[2],

        p_obrien = obrien_result$p_value,
        p_cv_tmle = cv_tmle_result$p_value,

        # CV-TMLE average weights and shrinkage diagnostics
        cv_alpha1 = cv_tmle_result$alpha_weights["alpha1"],
        cv_alpha2 = cv_tmle_result$alpha_weights["alpha2"],
        cv_lambda_mean = mean(cv_tmle_result$lambda_shrink, na.rm = TRUE),

        stringsAsFactors = FALSE
      )
    }

    do.call(rbind, out_list)
  }

  # ---- register parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl), add = TRUE)

  # Export helper functions to workers
  clusterExport(cl, c(
    "simulate_rank_change_trial",
    "tmle_multiplicity_test",
    "obrien_rank_sum_test",
    "cv_tmle_test",
    "mc_supremum_pvalue"
  ), envir = environment())

  # ---- run simulation via foreach/dopar
  final_results <- foreach(
    rep = 1:n_reps,
    .combine = rbind,
    .packages = c("tmle", "SuperLearner", "MASS", "dplyr"),
    .inorder = FALSE
  ) %dopar% {
    one_replication(
      rep = rep,
      scenarios = scenarios,
      n_total = n_total,
      alpha = alpha,
      covariate_vars = covariate_vars,
      tmle_libraries = tmle_libraries,
      cv_tmle_params = cv_tmle_params,
      seed = seed,
      tmle_verbose = tmle_verbose
    )
  }

  rownames(final_results) <- NULL

  if (verbose) {
    cat("\nSimulation completed.\n")
    cat("Total rows:", nrow(final_results), "\n")
  }

  final_results
}

# =============================================================================
# 4. RESULTS SUMMARY
# =============================================================================

summarize_four_method_results <- function(sim_results) {
  sim_results %>%
    group_by(scenario_name) %>%
    summarise(
      n_rows = n(),
      tmle_holm_rate = mean(tmle_holm_reject, na.rm = TRUE),
      tmle_hochberg_rate = mean(tmle_hochberg_reject, na.rm = TRUE),
      obrien_rate = mean(obrien_reject, na.rm = TRUE),
      cv_tmle_rate = mean(cv_tmle_reject, na.rm = TRUE),

      mean_corr = mean(observed_correlation, na.rm = TRUE),

      mean_cv_alpha1 = mean(cv_alpha1, na.rm = TRUE),
      mean_cv_alpha2 = mean(cv_alpha2, na.rm = TRUE),
      mean_cv_lambda = mean(cv_lambda_mean, na.rm = TRUE),

      .groups = "drop"
    ) %>%
    arrange(scenario_name)
}

# =============================================================================
# 5. MAIN EXECUTION
# =============================================================================

main_four_method_simulation <- function(
  n_reps = 1000,
  n_total = 60,
  seed = 202701,
  alpha = 0.025,
  save_results = FALSE,
  n_cores = parallel::detectCores() - 1
) {

  cat("=== FOUR-METHOD SIMULATION STUDY (REALISTIC DGP; CONTINUOUS CHANGES + RANK O'BRIEN) ===\n")
  cat("Methods: TMLE-Holm, TMLE-Hochberg, O'Brien(Unadj), Stabilized CV-TMLE\n")
  cat("Scenarios: global_null (distribution null) + table5_alt + weak_alt\n")
  cat("n_total:", n_total, " | n_reps:", n_reps, " | alpha:", alpha, "\n")
  cat("Using", n_cores, "cores for parallel execution\n\n")

  tmle_libs <- list(Q.SL.library = c("SL.glm"))

  start_time <- Sys.time()

  raw_results <- run_four_method_simulation(
    n_total = n_total,
    scenarios = list(
      list(name = "global_null", scenario = "null"),
      list(name = "table5_alt",  scenario = "alt"),
      list(name = "weak_alt",    scenario = "weak_alt")
    ),
    n_reps = n_reps,
    alpha = alpha,
    covariate_vars = c("age", "region", "sev_cat"),
    tmle_libraries = tmle_libs,

    cv_tmle_params = list(V = 10, C = 2, df_t = 36),

    seed = seed,
    verbose = TRUE,
    tmle_verbose = FALSE,
    n_cores = n_cores
  )

  end_time <- Sys.time()
  runtime <- end_time - start_time
  cat(sprintf("\nRuntime: %.2f minutes\n\n", as.numeric(runtime, units = "mins")))

  summary_results <- summarize_four_method_results(raw_results)

  cat("=== REJECTION RATES (Type I under global_null; Power under table5_alt and weak_alt) ===\n")
  print(summary_results %>%
          mutate(across(ends_with("_rate"), ~ round(.x, 3)),
                 mean_corr = round(mean_corr, 3),
                 mean_cv_alpha1 = round(mean_cv_alpha1, 3),
                 mean_cv_alpha2 = round(mean_cv_alpha2, 3),
                 mean_cv_lambda = round(mean_cv_lambda, 3)),
        n = Inf)

  if (save_results) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    write.csv(raw_results, paste0("realistic_sim_raw_", timestamp, ".csv"), row.names = FALSE)
    write.csv(summary_results, paste0("realistic_sim_summary_", timestamp, ".csv"), row.names = FALSE)
    cat("\nResults saved with timestamp:", timestamp, "\n")
  }

  list(
    raw_results = raw_results,
    summary_results = summary_results,
    runtime = runtime
  )
}


# =============================================================================
# 6. RUN SIMULATION
# =============================================================================
# NOTE: CV-TMLE is computationally expensive; start small, then scale up.
#
set.seed(202701)
res <- main_four_method_simulation(
  n_reps = 1000,
  n_total = 60,
  save_results = TRUE
)

str(res$raw_results)
print(res$summary_results)
