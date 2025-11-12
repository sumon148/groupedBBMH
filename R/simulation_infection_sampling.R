#' Simulate Infection Sampling with (Truncated) Beta-binomial adjusting Sensitivity and Minimum Positive Propensity
#'
#' Simulates contamination prevalence in batches using a (truncated) beta distribution,
#' generates counts of infected items per group, and models the effect of sample-based testing
#' under varying group sensitivities. Sampling within groups is done using hypergeometric draws,
#' reflecting real-world testing practices.
#'
#' This function helps compare the true contamination status in the population
#' versus what is observed through imperfect testing of sampled groups.
#'
#' @param alpha Numeric. Shape parameter \eqn{\alpha} of the beta distribution for prevalence.
#' @param beta Numeric. Shape parameter \eqn{\beta} of the beta distribution for prevalence.
#' @param cutoff Numeric (0 < cutoff ≤ 1). Lower bound for the truncated beta distribution.
#' @param D Integer. Number of targeted batches/replicates (i.e., number of simulated production units).
#' @param B Integer. Number of groups per batch.
#' @param M Integer. Number of items in each group.
#' @param m Integer. Number of items sampled from each group (for testing). If all of group items are tested, m=M.
#' @param b Integer. Number of groups sampled per batch.
#' @param sensitivity_values Numeric vector. A vector of sensitivity values (from 0 to 1)
#'        representing the probability that an infected group will be correctly identified.
#'
#' @return A named list of results, one element per sensitivity value. Each element is itself
#'         a list containing:
#' \describe{
#'   \item{`sensitivity`}{The sensitivity value used in that run.}
#'   \item{`T_Xi`}{True total number of contaminated items per batch.}
#'   \item{`T_Yi`}{True number of contaminated groups per batch (binary presence).}
#'   \item{`t_Xi`}{Total sampled contaminated items (with perfect sensitivity).}
#'   \item{`t_Yi`}{Number of sampled contaminated groups (with perfect sensitivity).}
#'   \item{`t_Xi_adjusted`}{Sampled contaminated items after applying sensitivity.}
#'   \item{`t_Yi_adjusted`}{Sampled contaminated groups after applying sensitivity.}
#'   \item{`comparison_tx`}{Comparison of true vs. detected contamination (items).}
#'   \item{`comparison_ty`}{Comparison of true vs. detected contamination (groups).}
#' }
#'
#' @details
#' The function simulates contamination as follows:
#' \enumerate{
#'   \item Draw a batch-level contamination prevalence from a (truncated) beta distribution.
#'   \item Simulate contaminated item counts in each group using a binomial distribution.
#'   \item Derive contaminated group status (0/1).
#'   \item Simulate sampling from groups using hypergeometric draws.
#'   \item Apply group-level sensitivity to reflect imperfect detection.
#'   \item Compare population vs. sample detection outcomes.
#' }
#'
#' Sensitivity values allow modeling of test performance degradation. The comparisons
#' (comparison_tx and comparison_ty) help quantify misclassification rates in sample-based testing.
#' Since the case study is based on the group-testing, comparison_ty is used for the case study.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' results <- simulate_infection_sampling(
#'   alpha = 0.005, beta = 4.5, cutoff = 0.01,
#'   D = 1000, B = 8000, M = 40, m = 5, b = 13,
#'   sensitivity_values = c(1, 0.8)
#' )
#'
#' # Inspect results for 90% sensitivity
#' str(results[["0.8"]])
#' }
#'
#' @seealso [rtrbeta()] for truncated beta distribution, [rhyper()] for hypergeometric sampling.
#'
#' @export
simulate_infection_sampling <- function(alpha, beta, cutoff, D, B, M, m, b, sensitivity_values) {

  #-----------------------------
  # Step 1: Generate batch-specific contamination prevalence
  #-----------------------------
  P_u <- rtrbeta(D, alpha, beta, cutoff)  # truncated beta distribution

  #-----------------------------
  # Step 2: Generate true number of infected cases per group (D x B)
  #-----------------------------
  X_ij <- matrix(NA, nrow = D, ncol = B)
  for (i in 1:D) {
    P_i <- P_u[i]
    X_ij[i, ] <- rbinom(B, size = M, prob = P_i)
  }

  # Total contaminated items per batch
  T_Xi <- apply(X_ij, 1, sum)

  #-----------------------------
  # Step 3: True infection status (binary) per group
  #-----------------------------
  Y_ij <- (X_ij > 0) * 1
  T_Yi <- apply(Y_ij, 1, sum)

  #-----------------------------
  # Step 4: Sampled infection status with perfect sensitivity
  #-----------------------------
  y_s_perfect <- matrix(NA, nrow = D, ncol = b)
  x_s_perfect <- matrix(NA, nrow = D, ncol = b)

  for (i in 1:D) {
    sampled_groups_i <- sample(c(1:B), b, replace = FALSE)

    if (m == M) {
      # Directly use sampled groups
      y_s_perfect[i, ] <- Y_ij[i, sampled_groups_i]
      x_s_perfect[i, ] <- X_ij[i, sampled_groups_i]
    } else {
      # Hypergeometric sampling from groups
      x_s_perfect[i, ] <- rhyper(b, X_ij[i, sampled_groups_i], M - X_ij[i, sampled_groups_i], m)
      y_s_perfect[i, ] <- (x_s_perfect[i, ] > 0) * 1
    }
  }

  # Summed sample counts
  t_Yi_perfect <- apply(y_s_perfect, 1, sum)
  t_Xi_perfect <- apply(x_s_perfect, 1, sum)

  #-----------------------------
  # Step 5: Apply imperfect sensitivity
  #-----------------------------
  results <- list()

  for (s in 1:length(sensitivity_values)) {
    sens <- sensitivity_values[s]

    y_s_adjusted <- matrix(NA, nrow = D, ncol = b)
    x_s_adjusted <- matrix(NA, nrow = D, ncol = b)

    for (i in 1:D) {
      # Sensitivity applied at the group level
      y_s_adjusted[i, ] <- rbinom(length(y_s_perfect[i, ]), size = y_s_perfect[i, ], prob = sens)

      # Keep counts only if group is detected positive
      x_s_adjusted[i, ] <- x_s_perfect[i, ] * I(y_s_adjusted[i, ] > 0)
    }

    # Summed sample counts (adjusted for sensitivity)
    t_Yi_adjusted <- apply(y_s_adjusted, 1, sum)
    t_Xi_adjusted <- apply(x_s_adjusted, 1, sum)

    #-----------------------------
    # Step 6: Compare population vs. sampled results
    #-----------------------------
    comparison_tx <- ifelse(T_Xi > 0 & t_Xi_adjusted > 0, 1,
                            ifelse(T_Xi == 0 & t_Xi_adjusted == 0, 2,
                                   ifelse(T_Xi > 0 & t_Xi_adjusted == 0, 3, 4)))

    comparison_ty <- ifelse(T_Yi > 0 & t_Yi_adjusted > 0, 1,
                            ifelse(T_Yi == 0 & t_Yi_adjusted == 0, 2,
                                   ifelse(T_Yi > 0 & t_Yi_adjusted == 0, 3, 4)))

    #-----------------------------
    # Step 7: Store results for current sensitivity
    #-----------------------------
    results[[as.character(sens)]] <- list(
      sensitivity    = sens,
      T_Xi           = T_Xi,
      t_Xi           = t_Xi_perfect,
      T_Yi           = T_Yi,
      t_Yi           = t_Yi_perfect,
      t_Yi_adjusted  = t_Yi_adjusted,
      t_Xi_adjusted  = t_Xi_adjusted,
      comparison_ty  = comparison_ty,
      comparison_tx  = comparison_tx
    )
  }

  return(results)
}



#' Run Repeated Infection Sampling Simulations
#'
#' Runs multiple simulations of infection prevalence and sampling using the
#' `simulate_infection_sampling()` function. For each sensitivity level, results are stored
#' across multiple simulated batches and iterations.
#'
#' @param alpha Numeric. Scalar or vector of shape parameters \eqn{\alpha} for the beta distribution.
#'        If scalar, it is repeated across simulations.
#' @param beta Numeric. Scalar or vector of shape parameters \eqn{\beta} for the beta distribution.
#'        If scalar, it is repeated across simulations.
#' @param cutoff Numeric. Scalar or vector of truncation points (0 < cutoff ≤ 1). If scalar, repeated.
#' @param D Integer. Number of batches per simulation.
#' @param B Integer. Number of groups per batch.
#' @param M Integer. Number of items in each group.
#' @param m Integer. Number of items sampled per group.
#' @param b Integer. Number of groups sampled per batch.
#' @param sensitivity_values Numeric vector. Vector of sensitivity values (between 0 and 1).
#' @param NSim Integer. Number of simulation runs.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @return A nested list containing simulation results for each sensitivity level. Each list includes:
#' \describe{
#'   \item{`Simulation_Ty`}{Matrix of true positive group counts (D x NSim).}
#'   \item{`Simulation_ty`}{Matrix of observed positive group counts after applying sensitivity.}
#'   \item{`Simulation_Tx`}{Matrix of true positive item counts per batch.}
#'   \item{`Simulation_tx`}{Matrix of observed positive item counts after sensitivity.}
#'   \item{`Simulation_Infection_Ty_Results`}{Comparison matrix for group-level detection.}
#'   \item{`Simulation_Infection_Tx_Results`}{Comparison matrix for item-level detection.}
#' }
#'
#' An additional element `Inputs` stores the function inputs for reference.
#'
#' @examples
#' \dontrun{
#' sim_results <- run_infection_simulations(
#'   alpha = 0.005, beta = 4.5, cutoff = 0.01,
#'   D = 1000, B = 8000, M = 40, m = 5, b = 13,
#'   sensitivity_values = c(1, 0.8),
#'   NSim = 100,
#'   seed = 123
#' )
#'
#' # View results for sensitivity = 0.8
#' str(sim_results[["0.8"]])
#' }
#'
#' @seealso [simulate_infection_sampling()]
#'
#' @export
run_infection_simulations <- function(alpha, beta, cutoff, D, B, M, m, b, sensitivity_values, NSim, seed) {

  # ---- FLOW OVERVIEW ----#
  # 1. Initialize → set seed, prepare result containers for each sensitivity value
  # 2. Expand inputs → ensure alpha, beta, cutoff have length = NSim
  # 3. Loop over simulations (i = 1...NSim):
  #      a. Call core simulation function (simulate_infection_sampling)
  #      b. Extract outputs for each sensitivity value
  #      c. Save results in corresponding matrices
  # 4. Attach input parameters to results
  # 5. Return results (nested list)
  #

  # ---- STEP 1: Initialize ----#
  set.seed(seed)                 # Reproducibility
  no.sim <- NSim                 # Number of simulations
  results <- list()              # Master results container

  # Create empty matrices for each sensitivity value
  for (s in sensitivity_values) {
    results[[as.character(s)]] <- list(
      Simulation_Ty = matrix(NA, D, no.sim),   # True infected groups (batch-level)
      Simulation_ty = matrix(NA, D, no.sim),   # Sampled infected groups
      Simulation_Tx = matrix(NA, D, no.sim),   # True infected units (item-level)
      Simulation_tx = matrix(NA, D, no.sim),   # Sampled infected units
      Simulation_Infection_Ty_Results = matrix(NA, D, no.sim), # Comparison outcomes (groups)
      Simulation_Infection_Tx_Results = matrix(NA, D, no.sim)  # Comparison outcomes (units)
    )
  }

  # ---- STEP 2: Expand inputs (if scalar) ----#
  alpha_vector  <- if (length(alpha)  == 1) rep(alpha,  no.sim) else alpha
  beta_vector   <- if (length(beta)   == 1) rep(beta,   no.sim) else beta
  cutoff_vector <- if (length(cutoff) == 1) rep(cutoff, no.sim) else cutoff

  # ---- STEP 3: Run simulation loop ----#
  for (i in 1:no.sim) {
    # (a) Call core simulation function → returns results for all sensitivities
    Sim <- simulate_infection_sampling(
      alpha = alpha_vector[i],
      beta = beta_vector[i],
      cutoff = cutoff_vector[i],
      D = D, B = B, M = M, m = m, b = b,
      sensitivity_values = sensitivity_values
    )

    # (b) Extract results for each sensitivity value
    for (s in sensitivity_values) {
      s_char <- as.character(s)
      results[[s_char]]$Simulation_Ty[, i] <- Sim[[s_char]]$T_Yi
      results[[s_char]]$Simulation_ty[, i] <- Sim[[s_char]]$t_Yi_adjusted
      results[[s_char]]$Simulation_Tx[, i] <- Sim[[s_char]]$T_Xi
      results[[s_char]]$Simulation_tx[, i] <- Sim[[s_char]]$t_Xi_adjusted
      results[[s_char]]$Simulation_Infection_Ty_Results[, i] <- Sim[[s_char]]$comparison_ty
      results[[s_char]]$Simulation_Infection_Tx_Results[, i] <- Sim[[s_char]]$comparison_tx
    }
  }

  # ---- STEP 4: Attach inputs ----#
  results$Inputs <- list(
    alpha = alpha,
    beta = beta,
    cutoff = cutoff,
    D = D,
    B = B,
    M = M,
    m = m,
    b = b,
    sensitivity_values = sensitivity_values,
    NSim = NSim
  )

  # ---- STEP 5: Return results ----#
  return(results)
}


#' Compute Expected Leakage and Probability of Leakage from Model-based Simulation Output
#'
#' This function computes batch-level and overall leakage statistics based on the output
#' of a simulation object generated by the package's core simulation functions
#' `simulate_infection_sampling` and `run_infection_simulations`.
#'
#' Leakage is defined as the number of undetected infections in a batch, particularly in
#' the presence of false negatives (FN). The function computes expected leakage, probability
#' of leakage, and pseudo-expected leakage (i.e., probability × median leakage).
#'
#' @param sim_object A simulation result object returned by the main simulation function,
#'   typically structured with elements: `Simulation_Tx`, `Simulation_tx`, and
#'   `Simulation_Infection_Ty_Results`.
#'
#' @details
#' The function calculates:
#' \itemize{
#'   \item \strong{E.Li}: Expected leakage per batch (mean of undetected infections)
#'   \item \strong{P.Li}: Probability of leakage (proportion of simulations with FN)
#'   \item \strong{Median.Li}: Median leakage, conditional on leakage > 0
#'   \item \strong{Pseudo.E.Li.median}: Product of P.Li and Median.Li
#'   \item \strong{Overall statistics}: Same metrics aggregated across all batches/simulations
#' }
#'
#' The function assumes false negatives are encoded as `3` in the `Simulation_Infection_Ty_Results`
#' matrix, and treats these as missed infections when computing leakage.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{E.Li}{Expected leakage per batch (mean undetected infections)}
#'   \item{P.Li}{Probability of leakage per batch}
#'   \item{Median.Li}{Median leakage per batch, conditional on leakage occurring}
#'   \item{Pseudo.E.Li.median}{Pseudo-expected leakage = P.Li × Median.Li}
#'   \item{Median.Li.Overall}{Overall median leakage across all batches/simulations}
#'   \item{P.Li.Overall}{Overall probability of leakage}
#'   \item{Pseudo.E.Li.Overall}{Overall pseudo-expected leakage}
#' }
#'
#'@seealso \code{\link{simulate_infection_sampling}}, \code{\link{run_infection_simulations}}
#'
#' @examples
#' # Example: Minimal simulation object with 5 batches and 5 simulations
#' Tx <- matrix(
#'   c(0, 13, 0, 8889, 0,
#'     0, 0, 0, 0, 0,
#'     0, 0, 0, 0, 5033,
#'     0, 0, 0, 0, 0,
#'     0, 0, 0, 0, 0),
#'   nrow = 5,
#'   ncol = 5,
#'   byrow = FALSE
#' )
#'
#' tx <- matrix(c(
#'   0, 0, 0, 0, 0,
#'   0, 0, 0, 0, 0,
#'   0, 0, 0, 0, 1,
#'   0, 0, 0, 0, 0,
#'   0, 0, 0, 0, 0
#' ), nrow = 5, ncol = 5, byrow = FALSE)
#'
#' FN <- matrix(c(
#'   2, 3, 2, 3, 2,
#'   2, 2, 2, 2, 2,
#'   2, 2, 2, 2, 1,
#'   2, 2, 2, 2, 2,
#'   2, 2, 2, 2, 2
#' ), nrow = 5, ncol = 5, byrow = FALSE)
#'
#' # Construct the simulation object
#' sim_obj <- list(
#'   Simulation_Tx = Tx,
#'   Simulation_tx = tx,
#'   Simulation_Infection_Ty_Results = FN
#' )
#'
#' # Compute leakage statistics
#' compute_simulated_leakage(sim_obj)
#' \dontrun{
#' sim <- run_infection_simulations(...)
#' leakage_stats <- compute_simulated_leakage(sim)
#' print(leakage_stats$E.Li)
#' }
#'
#' @export
compute_simulated_leakage <- function(sim_object) {
  # === 1. Extract relevant simulation matrices ===
  Tx  <- sim_object$Simulation_Tx                  # True number of infections
  tx  <- sim_object$Simulation_tx                  # Detected number of infections
  FN  <- sim_object$Simulation_Infection_Ty_Results == 3  # Indicator matrix for false negatives (FN == 3)

  # === 2. Calculate leakage: infections missed due to false negatives ===
  Tx_minus_tx <- Tx - tx                           # Undetected infections per simulation
  LKG <- Tx_minus_tx * FN                          # Apply only to false negatives (Expected Leakage: E(Li))

  # === 3. Compute expected leakage and leakage probability per batch ===
  E.Li <- apply(LKG, 2, mean)                      # Expected leakage per batch (average over simulations)
  P.Li <- apply(FN, 2, mean)                       # Probability of leakage per batch (FN rate)

  # === 4. Compute conditional leakage statistics (given leakage occurs) ===
  Median.Li <- apply(LKG, 2, function(x) median(x[x > 0]))  # Median leakage per batch when leakage > 0
  # Mean.Li <- apply(LKG, 1, function(x) mean(x[x > 0]))    # (Commented) Alternative: mean leakage when > 0

  # === 5. Compute pseudo-expected leakage (median * probability) ===
  Pseudo.E.Li.median <- P.Li * Median.Li           # Pseudo-expected leakage using median

  # === 6. Compute overall statistics (across all batches and simulations) ===
  LKG.vector <- as.vector(LKG)                     # Flatten the leakage matrix for overall stats
  Median.Li.Overall <- median(LKG.vector[LKG.vector > 0])  # Median leakage overall (only if leakage > 0)
  # Mean.Li.Overall <- mean(LKG.vector[LKG.vector > 0])     # (Commented) Mean leakage overall

  P.Li.Overall <- mean(FN)                         # Overall probability of leakage
  Pseudo.E.Li.Overall <- P.Li.Overall * Median.Li.Overall  # Overall pseudo-expected leakage
  # Pseudo.E.Li.mean.Overall <- P.Li.Overall * Mean.Li.Overall  # (Commented) Using mean instead
  E.Li.Overall <- mean(LKG)  # Mean of all values in the leakage matrix

  # === 7. Return results as a list ===
  results <- list(
    E.Li = E.Li,                                   # Expected leakage per batch
    P.Li = P.Li,                                   # Probability of leakage per batch
    Median.Li = Median.Li,                         # Median leakage per batch (given leakage)
    # Mean.Li = Mean.Li,                           # (Optional) Uncomment to return mean leakage
    Pseudo.E.Li.median = Pseudo.E.Li.median,       # Pseudo-expected leakage per batch (median-based)
    # Pseudo.E.Li.mean = Pseudo.E.Li.mean,         # (Optional) Uncomment to return mean-based
    Median.Li.Overall = Median.Li.Overall,         # Overall median leakage (given leakage)
    P.Li.Overall = P.Li.Overall,                   # Overall leakage probability
    E.Li.Overall = E.Li.Overall,                   # Overall expected leakage
    Pseudo.E.Li.Overall = Pseudo.E.Li.Overall      # Overall pseudo-expected leakage (median-based)
    # Pseudo.E.Li.mean.Overall = Pseudo.E.Li.mean.Overall  # (Optional)
  )

  return(results)
}





