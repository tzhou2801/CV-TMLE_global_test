# =============================================================================
# Simulation Study 1 (Section 6.1) for:
#   "A CV-TMLE global test approach to improve power in rare disease clinical
#    studies with multiple-component endpoints" (arXiv:2605.02851)
# =============================================================================
#
# This script implements a simulation study comparing FOUR
# global null hypothesis testing methods for two correlated outcomes.
#
# STABILIZATION P-VALUE:
# The stabilization p-value in the CV-TMLE method uses a Monte Carlo
# approximation to the null distribution of the supremum statistic
# T*_v = sup_{alpha} T_alpha (Section 4.5), which gives a properly
# calibrated p-value for the training-fold test.
#
# METHODS COMPARED:
# 1. Holm procedure (Note in this simple case with K=2 outcomes it is
#    equivalent to the Bonferroni correction.)
# 2. Hochberg procedure
# 3. O'Brien Global Test
# 4. Stabilized CV-TMLE with data-adaptive weights (Monte Carlo supremum p-value)
# Note methods 1-3 have also been applied to TMLE estimates (linear model for Q,
# known g0) for fair comparisons with CV-TMLE.
#
# Author: Tianyue Zhou, tianyue_zhou@berkeley.edu
# Date: March 2026
# =============================================================================

# Required libraries
library(dplyr)
library(tmle)
library(SuperLearner)
library(foreach)
library(doParallel)
library(MASS)  # for mvrnorm

# =============================================================================
# 1. OUTCOME GENERATION FUNCTION
# =============================================================================

genYs <- function(A, age, sex, 
                  betaA1, betaA2, 
                  beta.age = c(-0.1, -0.05), 
                  beta.sex = c(0.6, 0.3)) {
  n <- length(A)
  Y1 <- 1 + betaA1 * A + beta.age[1] * age + beta.sex[1] * sex + rnorm(n)
  Y2 <- 1 + betaA2 * A + beta.age[2] * age + beta.sex[2] * sex + rnorm(n)
  return(cbind(Y1 = Y1, Y2 = Y2))
}

# =============================================================================
# 2. SHARED TMLE COMPUTATION FOR HOLM, HOCHBERG, AND O'BRIEN
# =============================================================================

#' Compute TMLEs once and apply all three methods (Holm, Hochberg, O'Brien)
run_shared_tmle <- function(data, treatment_var, outcome_vars,
                            covariate_vars = c("age", "sex"),
                            alpha = 0.025,
                            one_sided = TRUE,
                            Q.SL.library = c("SL.glm")) {

  K <- length(outcome_vars)
  n <- nrow(data)
  A <- data[[treatment_var]]
  W <- data[, covariate_vars, drop = FALSE]

  # Known propensity score for RCT
  g0 <- rep(0.5, n)

  # Degrees of freedom
  df_holm_hoch <- n - 4
  df_obrien <- 30 # Small-sample correction for O'Brien from Logan and Tamhane

  # Run TMLE ONCE for each outcome
  psi <- numeric(K)
  se <- numeric(K)
  IC_matrix <- matrix(NA, nrow = n, ncol = K)

  for (k in 1:K) {
    tmle_k <- tmle(Y = data[[outcome_vars[k]]], A = A, W = W,
                   Q.SL.library = Q.SL.library,
                   g1W = g0,
                   family = "gaussian",
                   verbose = FALSE)

    psi[k] <- tmle_k$estimates$ATE$psi
    se[k] <- sqrt(tmle_k$estimates$ATE$var.psi)
    IC_matrix[, k] <- tmle_k$estimates$IC$IC.ATE
  }

  # Z-statistics
  z_stats <- psi / se

  # === Holm === 
  # In this simple case with K=2 outcomes, Holm is equivalent to Bonferroni
  if (one_sided) {
    p_values <- pt(z_stats, df = df_holm_hoch, lower.tail = FALSE)
  } else {
    p_values <- 2 * pt(abs(z_stats), df = df_holm_hoch, lower.tail = FALSE)
  }

  alpha_corrected <- alpha / K
  holm_reject <- any(p_values < alpha_corrected)

  # === Hochberg ===
  p_order <- order(p_values, decreasing = TRUE)
  p_sorted <- p_values[p_order]

  reject_individual <- logical(K)
  for (i in 1:K) {
    if (p_sorted[i] <= alpha / i) {
      for (j in i:K) {
        reject_individual[p_order[j]] <- TRUE
      }
      break
    }
  }
  hochberg_reject <- any(reject_individual)

  # === O'Brien ===
  R <- cor(IC_matrix)
  sum_z <- sum(z_stats)
  var_sum_z <- sum(R)
  t_obrien <- sum_z / sqrt(var_sum_z)
  p_obrien <- pt(t_obrien, df = df_obrien, lower.tail = FALSE)
  obrien_reject <- p_obrien < alpha

  return(list(
    estimates = psi,
    standard_errors = se,
    z_stats = z_stats,
    IC_matrix = IC_matrix,
    p_values = p_values,

    holm_reject = holm_reject,
    hochberg_reject = hochberg_reject,
    hochberg_reject_individual = reject_individual,

    obrien_reject = obrien_reject,
    obrien_pvalue = p_obrien,
    obrien_tstat = t_obrien,
    obrien_correlation = R[1, 2]
  ))
}

# =============================================================================
# 3. MONTE CARLO P-VALUE FOR SUPREMUM STATISTIC
# =============================================================================

#' Compute a Monte Carlo p-value for the supremum test statistic
#' T*_v = sup_{alpha in [0,1]} T_alpha, where
#' T_alpha = (alpha * Z_1 + (1-alpha) * Z_2) / sigma_alpha.
#'
#' Under H0, (Z_1, Z_2) ~ N(0, Sigma). We draw B samples from N(0, Sigma),
#' compute the supremum for each, and return the fraction exceeding observed_T.
#'
#' @param Sigma_mat 2x2 covariance matrix of influence curves (from training data)
#' @param observed_T The observed supremum statistic T*_v
#' @param B Number of Monte Carlo draws (default 5000)
#' @return Monte Carlo p-value P(T*_mc > observed_T)
mc_supremum_pvalue <- function(Sigma_mat, observed_T, B = 5000) {
  # Draw B samples from N(0, Sigma)
  Z_draws <- mvrnorm(n = B, mu = c(0, 0), Sigma = Sigma_mat)

  s11 <- Sigma_mat[1, 1]
  s12 <- Sigma_mat[1, 2]
  s22 <- Sigma_mat[2, 2]

  # For each draw, compute T* = sup_{alpha in [0,1]} T_alpha
  # T_alpha = (alpha * Z1 + (1-alpha) * Z2) / sqrt(alpha^2 * s11 + 2*alpha*(1-alpha)*s12 + (1-alpha)^2 * s22)
  # Evaluate at boundary values alpha=0, alpha=1
  T_at_0 <- Z_draws[, 2] / sqrt(s22)     # alpha = 0: weight entirely on Z2
  T_at_1 <- Z_draws[, 1] / sqrt(s11)     # alpha = 1: weight entirely on Z1

  # Find interior critical point for each draw using optimize()
  # The interior optimum of T(alpha) for given (Z1, Z2) can be found analytically
  # but optimize() is simpler and fast for 1D on B draws
  T_star_mc <- pmax(T_at_0, T_at_1)

  # Check interior: for each draw, optimize T(alpha) on (0,1)
  for (b in 1:B) {
    z1 <- Z_draws[b, 1]
    z2 <- Z_draws[b, 2]

    opt <- optimize(
      f = function(a) {
        num <- a * z1 + (1 - a) * z2
        denom <- sqrt(a^2 * s11 + 2 * a * (1 - a) * s12 + (1 - a)^2 * s22)
        return(-num / denom)  # minimize negative = maximize
      },
      interval = c(0, 1)
    )

    T_interior <- -opt$objective
    if (T_interior > T_star_mc[b]) {
      T_star_mc[b] <- T_interior
    }
  }

  # p-value: fraction of MC suprema exceeding observed value
  p_val <- mean(T_star_mc > observed_T)
  return(p_val)
}

# =============================================================================
# 4. CV-TMLE WITH DATA-ADAPTIVE WEIGHTS
# =============================================================================

# Stabilized CV-TMLE with gradient-based weight optimization and shrinkage
cv_tmle_test <- function(data, treatment_var, outcome1_var, outcome2_var,
                         covariate_vars = c("age", "sex"),
                         alpha = 0.025,
                         V = 10,  # number of folds
                         Q.SL.library = c("SL.glm"),
                         g.SL.library = c("SL.glm"),
                         verbose = FALSE,
                         seed = NULL) {

  # Degrees of freedom (small-sample correction for K=2 outcomes)
  df_cv <- 30
  
  n <- nrow(data)
  Y1 <- data[[outcome1_var]]
  Y2 <- data[[outcome2_var]]
  A <- data[[treatment_var]]
  W <- data[, covariate_vars, drop = FALSE]
  
  # Create folds
  folds <- sample(rep(1:V, length.out = n))
  
  # Storage for fold-specific results
  psi_star_v <- numeric(V)
  IC_all <- numeric(n)
  alpha_v_matrix <- matrix(NA, nrow = V, ncol = 2)
  
  # Process each fold
  for (v in 1:V) {
    if (verbose) cat("Processing fold", v, "of", V, "\n")
    
    # Split data
    train_idx <- which(folds != v)
    valid_idx <- which(folds == v)
    
    train_data <- data[train_idx, ]
    valid_data <- data[valid_idx, ]
    
    # Known propensity scores
    g0_train <- rep(0.5, nrow(train_data))
    g0_valid <- rep(0.5, nrow(valid_data))
    
    # Step 1: Compute individual TMLEs on training data
    tmle1_train <- tmle(Y = train_data[[outcome1_var]], 
                        A = train_data[[treatment_var]], 
                        W = train_data[, covariate_vars, drop = FALSE],
                        Q.SL.library = Q.SL.library,
                        g1W = g0_train,
                        family = "gaussian",
                        verbose = FALSE)
    
    tmle2_train <- tmle(Y = train_data[[outcome2_var]], 
                        A = train_data[[treatment_var]], 
                        W = train_data[, covariate_vars, drop = FALSE],
                        Q.SL.library = Q.SL.library,
                        g1W = g0_train,
                        family = "gaussian",
                        verbose = FALSE)
    
    # Extract point estimates
    psi1_star <- tmle1_train$estimates$ATE$psi
    psi2_star <- tmle2_train$estimates$ATE$psi
    
    # Step 2: Get influence curves and compute covariance
    IC1_train <- tmle1_train$estimates$IC$IC.ATE
    IC2_train <- tmle2_train$estimates$IC$IC.ATE
    
    sigma11 <- var(IC1_train)
    sigma22 <- var(IC2_train)
    sigma12 <- cov(IC1_train, IC2_train)
    
    # Step 3: L-BFGS-B optimization for optimal alpha
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
      grad_reduced <- dt_dalpha[1] - dt_dalpha[2]
      return(-grad_reduced)
    }

    optim_result <- tryCatch({
      optim(par = 0.5, fn = objective_fn, gr = gradient_fn,
            method = "L-BFGS-B", lower = 0, upper = 1,
            control = list(maxit = 200, factr = 1e7))
    }, error = function(e) {
      list(par = 0.5, value = -objective_fn(0.5))
    })

    alpha1_optimal <- optim_result$par
    alpha_optimal_temp <- c(max(0, alpha1_optimal), max(0, 1 - alpha1_optimal))
    alpha_optimal_temp <- alpha_optimal_temp / sum(alpha_optimal_temp)
    max_test_stat <- -optim_result$value

    # Step 3b: Shrinkage toward equal weights (stabilized CV-TMLE)
    alpha_ref <- c(0.5, 0.5)
    n_train_size <- n * (V - 1) / V
    max_test_stat_scaled <- sqrt(n_train_size) * max_test_stat
    # Monte Carlo p-value using null distribution of the supremum statistic
    p_n <- mc_supremum_pvalue(Sigma_mat, max_test_stat_scaled, B = 5000)
    C_tilde <- 0.25 * log(n)  # C=0.25, paper section 4.5
    shrinkage_weight <- min(C_tilde * p_n, 1)
    alpha_v <- (1 - shrinkage_weight) * alpha_optimal_temp +
      shrinkage_weight * alpha_ref

    alpha_v_matrix[v, ] <- alpha_v
    
    # Step 4: Apply optimal weights to validation data
    Y_bar_train <- alpha_v[1] * train_data[[outcome1_var]] + 
      alpha_v[2] * train_data[[outcome2_var]]
    Y_bar_valid <- alpha_v[1] * valid_data[[outcome1_var]] + 
      alpha_v[2] * valid_data[[outcome2_var]]
    
    # Fit initial outcome model for composite outcome
    train_model_df <- data.frame(
      Y_bar = Y_bar_train,
      A = train_data[[treatment_var]],
      train_data[, covariate_vars, drop = FALSE]
    )
    
    # Use linear model for simplicity and fair comparison, replace with SL in practice 
    lm_formula <- as.formula(paste("Y_bar ~ A +", paste(covariate_vars, collapse = " + ")))
    lm_composite <- lm(lm_formula, data = train_model_df)
    
    # Predict Q(A,W) for validation data
    valid_model_df <- data.frame(
      A = valid_data[[treatment_var]],
      valid_data[, covariate_vars, drop = FALSE]
    )
    
    Q_alpha_v <- predict(lm_composite, newdata = valid_model_df)
    
    # Predict Q(1,W) and Q(0,W)
    valid_model_df_1 <- valid_model_df
    valid_model_df_1$A <- 1
    Q_alpha_v_1 <- predict(lm_composite, newdata = valid_model_df_1)
    
    valid_model_df_0 <- valid_model_df
    valid_model_df_0$A <- 0
    Q_alpha_v_0 <- predict(lm_composite, newdata = valid_model_df_0)
    
    # Store for pooled targeting
    if (v == 1) {
      Y_bar_stacked <- Y_bar_valid
      Q_offset_stacked <- Q_alpha_v
      Q_1_stacked <- Q_alpha_v_1
      Q_0_stacked <- Q_alpha_v_0
      A_stacked <- valid_data[[treatment_var]]
      W_stacked <- valid_data[, covariate_vars, drop = FALSE]
      fold_id <- rep(v, length(valid_idx))
    } else {
      Y_bar_stacked <- c(Y_bar_stacked, Y_bar_valid)
      Q_offset_stacked <- c(Q_offset_stacked, Q_alpha_v)
      Q_1_stacked <- c(Q_1_stacked, Q_alpha_v_1)
      Q_0_stacked <- c(Q_0_stacked, Q_alpha_v_0)
      A_stacked <- c(A_stacked, valid_data[[treatment_var]])
      W_stacked <- rbind(W_stacked, valid_data[, covariate_vars, drop = FALSE])
      fold_id <- c(fold_id, rep(v, length(valid_idx)))
    }
  }
  
  # Step 5: Pooled targeting
  H_stacked <- (2 * A_stacked - 1) / 0.5  # Clever covariate for RCT with g0 = 0.5
  
  # Targeting step
  pooled_target <- lm(Y_bar_stacked ~ -1 + H_stacked, offset = Q_offset_stacked)
  epsilon <- coef(pooled_target)
  
  # Update predictions
  Q_star_stacked <- Q_offset_stacked + epsilon * H_stacked
  Q_star_1_stacked <- Q_1_stacked + epsilon * 2      # H=2 when A=1
  Q_star_0_stacked <- Q_0_stacked + epsilon * (-2)    # H=-2 when A=0
  
  # Step 6: Compute fold-specific TMLEs and influence curves
  for (v in 1:V) {
    fold_mask <- (fold_id == v)
    
    # Get targeted predictions for this fold
    Q_star_1_fold <- Q_star_1_stacked[fold_mask]
    Q_star_0_fold <- Q_star_0_stacked[fold_mask]
    Q_star_fold <- Q_star_stacked[fold_mask]
    Y_bar_fold <- Y_bar_stacked[fold_mask]
    H_fold <- H_stacked[fold_mask]
    
    # TMLE estimate for fold v
    psi_star_v[v] <- mean(Q_star_1_fold - Q_star_0_fold)
    
    # Influence curves
    IC_fold <- H_fold * (Y_bar_fold - Q_star_fold) + 
      (Q_star_1_fold - Q_star_0_fold) - psi_star_v[v]
    
    # Store in correct positions
    IC_all[which(folds == v)] <- IC_fold
  }
  
  # Step 7: Final CV-TMLE estimate and inference
  psi_cv_tmle <- mean(psi_star_v)
  sigma2_cv <- var(IC_all)
  sigma_cv <- sqrt(sigma2_cv)
  
  # Calculate average alpha across all folds
  alpha_avg <- colMeans(alpha_v_matrix)
  names(alpha_avg) <- c("alpha1", "alpha2")
  
  # Test statistic and decision
  t_stat <- sqrt(n) * psi_cv_tmle / sigma_cv
  p_value <- pt(t_stat, df = df_cv, lower.tail = FALSE)  # One-sided test

  # Return results
  return(list(
    reject_h0 = t_stat > qt(1 - alpha, df = df_cv),
    p_value = p_value,
    psi_cv_tmle = psi_cv_tmle,
    sigma_cv = sigma_cv,
    t_stat = t_stat,
    alpha_weights = alpha_avg,
    alpha_weights_matrix = alpha_v_matrix
  ))
}

# =============================================================================
# 5. COMPREHENSIVE SIMULATION PIPELINE
# =============================================================================

run_four_method_simulation <- function(n_total = 50,
                                       scenarios = list(
                                         list(betaA1 = 0.0, betaA2 = 0.0, name = "Y00_type_1_error"),
                                         list(betaA1 = 1.0, betaA2 = 0.0, name = "Y10_strong_Y1_only"),
                                         list(betaA1 = 0.5, betaA2 = 0.5, name = "Y55_equal_moderate"),
                                         list(betaA1 = 0.8, betaA2 = 0.2, name = "Y82_strong_Y1_weak_Y2")
                                       ),
                                       age_effects = list(
                                         list(effects = c(-0.1, -0.05), name = "moderate"),
                                         list(effects = c(-0.5, -0.2), name = "strong")
                                       ),
                                       n_reps = 500,
                                       alpha = 0.025,
                                       covariate_vars = c("age", "sex"),
                                       tmle_libraries = list(
                                         Q.SL.library = c("SL.glm"),
                                         g.SL.library = c("SL.glm")
                                       ),
                                       seed = 202701,
                                       n_cores = parallel::detectCores() - 1,
                                       verbose = TRUE,
                                       tmle_verbose = FALSE) {

  total_combinations <- length(scenarios) * length(age_effects) * n_reps

  if (verbose) {
    cat("Starting FOUR-METHOD simulation study...\n")
    cat("Methods: TMLE-Holm, TMLE-Hochberg, TMLE-O'Brien, CV-TMLE\n")
    cat("Total combinations:", total_combinations, "\n")
    cat("Using", n_cores, "cores for parallel execution\n\n")
  }

  # Register parallel backend
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  on.exit(stopCluster(cl))

  # Export helper functions to workers
  clusterExport(cl, c("genYs", "run_shared_tmle", "cv_tmle_test",
                       "mc_supremum_pvalue"),
                envir = environment(genYs))

  final_results <- foreach(
    rep = 1:n_reps,
    .combine = rbind,
    .packages = c("tmle", "SuperLearner", "dplyr", "MASS"),
    .inorder = FALSE
  ) %dopar% {
    # Per-replication seed for reproducibility
    set.seed(seed + rep)

    # Generate baseline data
    age <- sample(5:18, n_total, replace = TRUE)
    sex <- rbinom(n_total, 1, plogis(0.3))
    A <- rbinom(n_total, 1, 0.5)

    rep_results <- vector("list", length(scenarios) * length(age_effects))
    idx <- 1

    for (s_idx in seq_along(scenarios)) {
      scenario <- scenarios[[s_idx]]

      for (a_idx in seq_along(age_effects)) {
        age_effect <- age_effects[[a_idx]]

        # Generate outcomes
        outcomes <- genYs(A, age, sex,
                          betaA1 = scenario$betaA1,
                          betaA2 = scenario$betaA2,
                          beta.age = age_effect$effects)

        sim_data <- data.frame(A = A, age = age, sex = sex,
                               Y1 = outcomes[, 1], Y2 = outcomes[, 2])

        # Run shared TMLE once for Holm, Hochberg, O'Brien
        shared_result <- run_shared_tmle(sim_data, "A", c("Y1", "Y2"),
                                         covariate_vars = covariate_vars,
                                         alpha = alpha,
                                         Q.SL.library = tmle_libraries$Q.SL.library)

        # Run CV-TMLE (separate, uses cross-validation)
        cv_tmle_result <- cv_tmle_test(sim_data, "A", "Y1", "Y2",
                                       covariate_vars = covariate_vars,
                                       alpha = alpha,
                                       V = 10,
                                       Q.SL.library = tmle_libraries$Q.SL.library,
                                       g.SL.library = tmle_libraries$g.SL.library,
                                       verbose = FALSE)

        # Store results
        rep_results[[idx]] <- data.frame(
          replication = rep,
          scenario_name = scenario$name,
          age_effect_name = age_effect$name,
          betaA1 = scenario$betaA1,
          betaA2 = scenario$betaA2,
          beta_age1 = age_effect$effects[1],
          beta_age2 = age_effect$effects[2],
          observed_correlation = cor(sim_data$Y1, sim_data$Y2),

          # Method rejections
          tmle_holm_reject = shared_result$holm_reject,
          tmle_hochberg_reject = shared_result$hochberg_reject,
          tmle_obrien_reject = shared_result$obrien_reject,
          cv_tmle_reject = cv_tmle_result$reject_h0,

          # P-values
          p_value_Y1 = shared_result$p_values[1],
          p_value_Y2 = shared_result$p_values[2],
          p_value_obrien = shared_result$obrien_pvalue,
          p_value_cv_tmle = cv_tmle_result$p_value,

          # Estimates (from shared TMLE, computed once)
          estimate_Y1 = shared_result$estimates[1],
          estimate_Y2 = shared_result$estimates[2],

          # CV-TMLE alpha weights
          cv_tmle_alpha1 = cv_tmle_result$alpha_weights["alpha1"],
          cv_tmle_alpha2 = cv_tmle_result$alpha_weights["alpha2"],

          # O'Brien correlation
          obrien_correlation = shared_result$obrien_correlation,

          stringsAsFactors = FALSE
        )

        idx <- idx + 1
      }
    }

    do.call(rbind, rep_results)
  }

  rownames(final_results) <- NULL

  if (verbose) {
    cat("\nFour-method simulation completed!\n")
    cat("Total rows:", nrow(final_results), "\n")
  }

  return(final_results)
}

# =============================================================================
# 6. RESULTS SUMMARY FUNCTION
# =============================================================================

summarize_four_method_results <- function(sim_results, alpha = 0.025) {
  summary_stats <- sim_results %>%
    group_by(scenario_name, age_effect_name, betaA1, betaA2, beta_age1, beta_age2) %>%
    summarise(
      n_replications = n(),
      
      # Power calculations for all four methods
      tmle_holm_power = mean(tmle_holm_reject, na.rm = TRUE),
      tmle_hochberg_power = mean(tmle_hochberg_reject, na.rm = TRUE),
      tmle_obrien_power = mean(tmle_obrien_reject, na.rm = TRUE),
      cv_tmle_power = mean(cv_tmle_reject, na.rm = TRUE),

      # Average CV-TMLE alpha weights
      mean_cv_tmle_alpha1 = mean(cv_tmle_alpha1, na.rm = TRUE),
      mean_cv_tmle_alpha2 = mean(cv_tmle_alpha2, na.rm = TRUE),
      sd_cv_tmle_alpha1 = sd(cv_tmle_alpha1, na.rm = TRUE),
      sd_cv_tmle_alpha2 = sd(cv_tmle_alpha2, na.rm = TRUE),
      
      # Correlation characteristics
      mean_observed_correlation = mean(observed_correlation, na.rm = TRUE),
      
      .groups = 'drop'
    ) %>%
    arrange(scenario_name, age_effect_name)
  
  return(summary_stats)
}

# =============================================================================
# 7. MAIN EXECUTION FUNCTION
# =============================================================================

main_four_method_simulation <- function(n_reps = 500, seed = 202701,
                                        n_cores = parallel::detectCores() - 1,
                                        save_results = FALSE) {
  
  cat("=== FOUR-METHOD SIMULATION STUDY (MC Supremum p-value) ===\n")
  cat("Comparing: TMLE-Holm, TMLE-Hochberg, TMLE-O'Brien, CV-TMLE\n")
  cat("CV-TMLE uses Monte Carlo calibrated p-value for supremum statistic\n\n")
  
  # Define scenarios
  treatment_scenarios <- list(
    list(betaA1 = 0.0, betaA2 = 0.0, name = "Y00_type_1_error"),
    list(betaA1 = 1.0, betaA2 = 0.0, name = "Y10_strong_Y1_only"),
    list(betaA1 = 0.5, betaA2 = 0.5, name = "Y55_equal_moderate"),
    list(betaA1 = 0.8, betaA2 = 0.2, name = "Y82_strong_Y1_weak_Y2")
  )

  age_effect_scenarios <- list(
    list(effects = c(-0.1, -0.05), name = "moderate_age_effects")#,
    #list(effects = c(-0.5, -0.2), name = "strong_age_effects")
  )
  
  tmle_libs <- list(Q.SL.library = c("SL.glm"), g.SL.library = c("SL.glm"))
  
  cat("Running simulation with", n_reps, "replications...\n")
  
  start_time <- Sys.time()
  
  raw_results <- run_four_method_simulation(
    n_total = 50, # Sample size
    scenarios = treatment_scenarios,
    age_effects = age_effect_scenarios,
    n_reps = n_reps,
    alpha = 0.025,
    covariate_vars = c("age", "sex"),
    tmle_libraries = tmle_libs,
    seed = seed,
    n_cores = n_cores,
    verbose = TRUE,
    tmle_verbose = FALSE
  )
  
  end_time <- Sys.time()
  runtime <- end_time - start_time
  
  cat(sprintf("\nRuntime: %.2f minutes\n\n", as.numeric(runtime, units = "mins")))
  
  # Generate summary
  summary_results <- summarize_four_method_results(raw_results)
  
  # Display results
  cat("=== POWER COMPARISON RESULTS ===\n")
  power_table <- summary_results %>%
    dplyr::select(scenario_name, age_effect_name, 
           tmle_holm_power, tmle_hochberg_power, 
           tmle_obrien_power, cv_tmle_power,
           mean_observed_correlation) %>%
    mutate(across(ends_with("power"), ~ round(.x, 3)),
           mean_observed_correlation = round(mean_observed_correlation, 3))
  
  print(power_table, n = Inf)
  
  # Display Type I error rates (Y00 scenario)
  cat("\n=== TYPE I ERROR RATES (Y00 Scenario) ===\n")
  type1_table <- summary_results %>%
    filter(scenario_name == "Y00_type_1_error") %>%
    dplyr::select(age_effect_name, 
           tmle_holm_power, tmle_hochberg_power, 
           tmle_obrien_power, cv_tmle_power) %>%
    rename(
      tmle_holm_type1 = tmle_holm_power,
      tmle_hochberg_type1 = tmle_hochberg_power,
      tmle_obrien_type1 = tmle_obrien_power,
      cv_tmle_type1 = cv_tmle_power
    ) %>%
    mutate(across(ends_with("type1"), ~ round(.x, 3)))
  
  print(type1_table)
  
  # Display CV-TMLE alpha weights
  cat("\n=== CV-TMLE ALPHA WEIGHT STATISTICS ===\n")
  alpha_stats_table <- summary_results %>%
    dplyr::select(scenario_name, age_effect_name,
           mean_cv_tmle_alpha1, mean_cv_tmle_alpha2,
           sd_cv_tmle_alpha1, sd_cv_tmle_alpha2) %>%
    mutate(across(starts_with("mean_"), ~ round(.x, 3)),
           across(starts_with("sd_"), ~ round(.x, 3)))

  print(alpha_stats_table, n = Inf)
  
  # Save results if requested
  if (save_results) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    write.csv(raw_results, paste0("simple_sim_supremum_raw_", timestamp, ".csv"), row.names = FALSE)
    write.csv(summary_results, paste0("simple_sim_supremum_summary_", timestamp, ".csv"), row.names = FALSE)
    cat("\nResults saved with timestamp:", timestamp, "\n")
  }
  
  return(list(
    raw_results = raw_results,
    summary_results = summary_results,
    runtime = runtime
  ))
}

# =============================================================================
# 8. RUN SIMULATION
# =============================================================================

# Execute the simulation
cat("Starting four-method simulation...\n")
cat("This will compare: Holm, Hochberg, O'Brien, and CV-TMLE\n\n")

results <- main_four_method_simulation(
  n_reps = 1000,      # Number of replications
  seed = 202701,   # Random seed for reproducibility
  n_cores = 7,
  save_results = TRUE  # Set to TRUE to save CSV files
)

cat("\n=== SIMULATION COMPLETED ===\n")
cat("Access results via:\n")
cat("- results$summary_results (power comparison and statistics)\n")
cat("- results$raw_results (detailed simulation data)\n\n")