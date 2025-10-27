#' Dunn-Smyth Residual Diagnostics for Beta-Binomial Models
#'
#' Computes Dunn-Smyth randomized quantile residuals for beta-binomial models,
#' including approximate and simulation-based versions depending on whether
#' posterior samples for \code{alpha} and \code{beta} are vectors or scalars.
#' Useful for diagnosing model fit in pooled testing or overdispersed count models.
#'
#' @param ObservedResponse A numeric vector of observed counts (length \eqn{D}).
#' @param size Integer. The number of trials in the binomial distribution (e.g., pool size).
#' @param group_size Integer. The size of the pooling group (used to adjust \code{beta}).
#' @param alpha A scalar or vector of \eqn{\alpha} values. If scalar, assumed fixed; if vector, assumed posterior samples.
#' @param beta A scalar or vector of \eqn{\beta} values. Same behavior as \code{alpha}.
#' @param approximate.model Logical. If \code{TRUE} (default), uses an approximate model where
#'        \eqn{\Phi_i \sim \text{Beta}(\alpha, \beta / \text{group\_size})}. If \code{FALSE}, simulates the pooled response using adjusted probabilities.
#'
#' @return A named list containing:
#' \describe{
#'   \item{\code{DS}}{Dunn-Smyth residuals on \[0,1\] scale.}
#'   \item{\code{normalized_DS}}{Dunn-Smyth residuals transformed via \code{qnorm()}.}
#'   \item{\code{KS}}{The result of a Kolmogorov-Smirnov test comparing residuals to Uniform(0,1).}
#'   \item{\code{DS_hist}}{A ggplot2 histogram of the residuals.}
#'   \item{\code{DS_QQ_norm}}{A normal Q-Q plot of the transformed residuals, stratified by observed counts.}
#' }
#'
#' @importFrom VGAM pbetabinom.ab
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_histogram labs annotate theme_minimal theme scale_color_manual
#' @importFrom stats ks.test pbinom runif qqnorm qnorm
#' @examples
#' \dontrun{
#'   set.seed(123)
#'   obs <- rbinom(100, size = 10, prob = 0.3)
#'   result <- DunnSmythTestBB(obs, size = 10, group_size = 5, alpha = 2, beta = 5)
#'   print(result$DS_QQ_norm)
#'   print(result$DS_hist)
#' }
#' @export
DunnSmythTestBB <- function(ObservedResponse, size, group_size, alpha, beta, approximate.model = TRUE) {

  ObservedResponse <- sort(ObservedResponse)
  D <- length(ObservedResponse)

  # Handle scalar or vector alpha/beta
  alpha_hat <- if (length(alpha) == 1) alpha else as.vector(alpha)
  beta_hat  <- if (length(beta) == 1) beta  else as.vector(beta)

  if (length(alpha) == 1 || length(beta) == 1) {
    # Scalar alpha/beta: Compute CDF directly
    size_vector <- rep(size, D)
    alpha_vec  <- rep(alpha_hat, D)
    beta_vec   <- rep(beta_hat, D)
    valid_idx  <- ObservedResponse > 0

    F_ymin1 <- rep(0, D)

    if (approximate.model) {
      F_ymin1[valid_idx] <- pbetabinom.ab(ObservedResponse[valid_idx] - 1, size_vector[valid_idx], alpha_vec[valid_idx], beta_vec[valid_idx] / group_size)
      F_y <- pbetabinom.ab(ObservedResponse, size_vector, alpha_vec, beta_vec / group_size)
    } else {
      N_samples <- 1000
      probs <- rowMeans(matrix(rbeta(D * N_samples, unique(alpha_vec), unique(beta_vec)), ncol = N_samples))
      probs_m <- 1 - (1 - probs)^group_size
      F_ymin1[valid_idx] <- pbinom(ObservedResponse[valid_idx] - 1, size_vector[valid_idx], probs_m[valid_idx])
      F_y <- pbinom(ObservedResponse, size_vector, probs_m)
    }
  } else {
    # Posterior samples: Simulate CDF via empirical quantiles
    S <- length(alpha_hat)
    sim_y <- matrix(nrow = D, ncol = S)

    for (s in seq_len(S)) {
      if (approximate.model) {
        probs <- rbeta(D, alpha_hat[s], beta_hat[s] / group_size)
        sim_y[, s] <- rbinom(D, size = size, prob = probs)
      } else {
        probs <- rbeta(D, alpha_hat[s], beta_hat[s])
        probs_m <- 1 - (1 - probs)^group_size
        sim_y[, s] <- rbinom(D, size = size, prob = probs_m)
      }
    }

    F_ymin1 <- rowMeans(sim_y < ObservedResponse)
    F_y <- rowMeans(sim_y <= ObservedResponse)
  }

  # Compute Dunn-Smyth residuals
  epsilon <- 1e-12
  F_ymin1 <- pmin(pmax(F_ymin1, epsilon), 1 - epsilon)
  F_y <- pmin(pmax(F_y, epsilon), 1 - epsilon)

  DS <- runif(D, F_ymin1, F_y)
  normalized_DS <- qnorm(pmin(pmax(DS, 1e-10), 1 - 1e-10))
  KS <- ks.test(DS, "punif")

  # Normal Q-Q plot
  qq_data <- data.frame(
    residuals = normalized_DS,
    group = ifelse(ObservedResponse > 0, "Non-zero", "Zero")
  )

  qq <- qqnorm(qq_data$residuals, plot.it = FALSE)

  qq_df <- data.frame(
    theoretical = qq$x,
    sample = qq$y,
    group = qq_data$group
  )

  DS_QQ_norm <- ggplot(qq_df, aes(x = theoretical, y = sample, color = group)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1.2) +
    labs(
      title = "Normal Q-Q Plot: Dunn-Smyth scaled residuals",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles",
      color = "Group"
    ) +
    scale_color_manual(values = c("Zero" = "grey", "Non-zero" = "blue")) +
    annotate("text", x = min(qq_df$theoretical), y = max(qq_df$sample),
             label = paste("KS test for Uniform: p =", round(KS$p.value, 4)),
             hjust = 0, vjust = 1, size = 4) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )

  # Histogram of DS residuals
  breaks <- seq(0, 1, length.out = 31)
  DS_hist <- ggplot(data.frame(residuals = DS), aes(x = residuals)) +
    geom_histogram(breaks = breaks, fill = "lightblue", color = "black", alpha = 0.7) +
    labs(
      title = "Distribution of Dunn-Smyth Residuals",
      x = "Residuals",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )

  return(list(
    DS = DS,
    normalized_DS = normalized_DS,
    KS = KS,
    DS_hist = DS_hist,
    DS_QQ_norm = DS_QQ_norm
  ))
}


#' Dunn-Smyth Residuals for Beta-Binomial Models (Custom Plotting)
#'
#' Computes Dunn-Smyth randomized quantile residuals for a beta-binomial model.
#' This custom version produces simpler diagnostic plots without titles or KS annotations.
#'
#' @param ObservedResponse A numeric vector of observed responses.
#' @param size Integer. Number of binomial trials per unit (e.g., pool size).
#' @param group_size Integer. Size of pooling group. Used to adjust probability calculation.
#' @param alpha Scalar or vector of alpha values. If vector, assumed to be posterior samples.
#' @param beta Scalar or vector of beta values. If vector, assumed to be posterior samples.
#' @param approximate.model Logical. If \code{TRUE} (default), uses approximate beta-binomial model with \eqn{\beta / \text{group\_size}}. If \code{FALSE}, simulates pooling process.
#'
#' @return A list containing:
#' \describe{
#'   \item{DS}{Raw Dunn-Smyth residuals on \[0,1\] scale.}
#'   \item{normalized_DS}{Residuals transformed to the standard normal scale.}
#'   \item{KS}{Kolmogorov-Smirnov test comparing residuals to Uniform(0,1).}
#'   \item{DS_hist}{A \code{ggplot2} histogram of residuals (no title).}
#'   \item{DS_QQ_norm}{Q-Q plot of normal Dunn-Smyth residuals, stratified by response group (no title, no annotation).}
#' }
#'
#' @importFrom VGAM pbetabinom.ab
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_histogram labs annotate theme_minimal theme scale_color_manual
#' @importFrom stats ks.test pbinom runif qqnorm qnorm
#'
#' @examples
#' \dontrun{
#' obs <- rbinom(100, size = 10, prob = 0.25)
#' result <- DunnSmythTestBB_Custom(obs, size = 10, group_size = 5, alpha = 2, beta = 4)
#' result$DS_QQ_norm
#' result$DS_hist
#' }
#' @export
DunnSmythTestBB_Custom <- function(ObservedResponse, size, group_size, alpha, beta, approximate.model = TRUE) {

  ObservedResponse <- sort(ObservedResponse)
  D <- length(ObservedResponse)

  alpha_hat <- if (length(alpha) == 1) alpha else as.vector(alpha)
  beta_hat <- if (length(beta) == 1) beta else as.vector(beta)

  if (length(alpha) == 1 || length(beta) == 1) {
    size_vector <- rep(size, D)
    alpha_vec <- rep(alpha_hat, D)
    beta_vec <- rep(beta_hat, D)
    valid_idx <- ObservedResponse > 0
    F_ymin1 <- rep(0, D)

    if (approximate.model) {
      F_ymin1[valid_idx] <- pbetabinom.ab(ObservedResponse[valid_idx] - 1, size_vector[valid_idx], alpha_vec[valid_idx], beta_vec[valid_idx] / group_size)
      F_y <- pbetabinom.ab(ObservedResponse, size_vector, alpha_vec, beta_vec / group_size)
    } else {
      N_samples <- 1000
      probs <- rowMeans(matrix(rbeta(D * N_samples, unique(alpha_vec), unique(beta_vec)), ncol = N_samples))
      probs_m <- 1 - (1 - probs)^group_size
      F_ymin1[valid_idx] <- pbinom(ObservedResponse[valid_idx] - 1, size_vector[valid_idx], probs_m[valid_idx])
      F_y <- pbinom(ObservedResponse, size_vector, probs_m)
    }

  } else {
    S <- length(alpha_hat)
    sim_y <- matrix(nrow = D, ncol = S)

    for (s in seq_len(S)) {
      if (approximate.model) {
        probs <- rbeta(D, alpha_hat[s], beta_hat[s] / group_size)
        sim_y[, s] <- rbinom(D, size = size, prob = probs)
      } else {
        probs <- rbeta(D, alpha_hat[s], beta_hat[s])
        probs_m <- 1 - (1 - probs)^group_size
        sim_y[, s] <- rbinom(D, size = size, prob = probs_m)
      }
    }

    F_ymin1 <- rowMeans(sim_y < ObservedResponse)
    F_y <- rowMeans(sim_y <= ObservedResponse)
  }

  epsilon <- 1e-12
  F_ymin1 <- pmin(pmax(F_ymin1, epsilon), 1 - epsilon)
  F_y <- pmin(pmax(F_y, epsilon), 1 - epsilon)

  DS <- runif(D, F_ymin1, F_y)
  normalized_DS <- qnorm(pmin(pmax(DS, 1e-10), 1 - 1e-10))
  KS <- ks.test(DS, "punif")

  qq_data <- data.frame(
    residuals = normalized_DS,
    group = ifelse(ObservedResponse > 0, "Non-zero", "Zero")
  )

  qq <- qqnorm(qq_data$residuals, plot.it = FALSE)

  qq_df <- data.frame(
    theoretical = qq$x,
    sample = qq$y,
    group = qq_data$group
  )

  DS_QQ_norm <- ggplot(qq_df, aes(x = theoretical, y = sample, color = group)) +
    geom_point(alpha = 0.7) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1.2) +
    labs(
      x = "Theoretical Quantiles",
      y = "Sample Quantiles",
      color = "Group"
    ) +
    scale_color_manual(values = c("Zero" = "grey", "Non-zero" = "blue")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 8)
    )

  DS_hist <- ggplot(data.frame(residuals = DS), aes(x = residuals)) +
    geom_histogram(breaks = seq(0, 1, length.out = 31), fill = "lightblue", color = "black", alpha = 0.7) +
    labs(
      x = "Residuals",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )

  list(
    DS = DS,
    normalized_DS = normalized_DS,
    KS = KS,
    DS_hist = DS_hist,
    DS_QQ_norm = DS_QQ_norm
  )
}
