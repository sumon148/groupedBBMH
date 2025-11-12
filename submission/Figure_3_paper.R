# --------------------------------------------------------------------------- #
# Load Required Libraries
# --------------------------------------------------------------------------- #
library(groupedBBMH)  # Custom package for grouped beta-binomial MH modeling
library(dplyr)        # Data manipulation
library(xtable)       # LaTeX table generation
library(bayesplot)    # Diagnostic and visualization tools for MCMC
library(coda)         # MCMC diagnostics and summaries
library(ggplot2)      # General plotting
library(ggpubr)       # Publication-ready plots and arrangement tools
library(MASS)         # For multivariate normal distribution functions
library(VGAM)         # Functions for vector generalized linear and additive models; used here for beta-binomial distribution functions (e.g., pbetabinom.ab)




#-------------------------------------------------------------@
# Figure 3: JASA PAper
# Dunn-Smyth Residual Goodness-of-Fit: MH Algorithm
# Using Seafood data
# Exact BB Model with \Delta=0.80: MH Algorithm
# Exact BB Model with \Delta ~ U(0.75, 0.85): MH Algorithm
#-------------------------------------------------------------@



DunnSmythTestBB <- function(ObservedResponse,size,group_size,alpha,beta,approximate.model=TRUE) {
  # ObservedResponse - Vector of dimension D
  # SimulatedResponse - Matrix of Dimension D X S
  # D for the length of observations
  # S for the number of simulated response
  # group_size for pool size
  # approximate.model for approximate distribution of Phi_i ~ Beta (alpha,beta/group_size)


  library(statmod)  # For qnorm and randomization
  library(VGAM)     # For pbetabinom.ab function
  library(ggplot2)

  ObservedResponse <- sort(ObservedResponse)

  D <- length(ObservedResponse) # Number of observed points
  # S <- dim(SimulatedResponse)[2]
  if (length(alpha)==1) {
    alpha_hat = alpha
  } else {alpha_hat=as.vector(alpha)}

  if (length(beta)==1) {
    beta_hat = beta
  } else {beta_hat=as.vector(beta)}

  if (length(alpha)==1 || length(beta)==1) {
    size_vector <- rep(size,D)
    alpha_hat_vector <- rep(alpha_hat,D)
    beta_hat_vector <- rep(beta_hat,D)
    valid_idx <- ObservedResponse > 0  # Identify indices where y_obs > 0
    F_ymin1 <- rep(0,length(ObservedResponse))

    if (approximate.model) {
      F_ymin1[valid_idx] <- pbetabinom.ab(ObservedResponse[valid_idx] - 1, size_vector[valid_idx], alpha_hat_vector[valid_idx], beta_hat_vector[valid_idx]/group_size)
      F_y <- pbetabinom.ab(ObservedResponse, size_vector, alpha_hat_vector, beta_hat_vector/group_size)
    } else {

      # probs <- rbeta(D,unique(alpha_hat_vector), unique(beta_hat_vector))
      # probs_m <- 1 - (1-probs)^group_size

      N_samples <- 1000
      probs <- rowMeans(matrix(rbeta(D * N_samples, unique(alpha_hat_vector), unique(beta_hat_vector)), ncol = N_samples))
      probs_m <- 1 - (1 - probs)^group_size

      F_ymin1[valid_idx] <- pbinom(ObservedResponse[valid_idx] - 1, size_vector[valid_idx], probs_m[valid_idx])
      F_y <- pbinom(ObservedResponse, size_vector, probs_m)
    }
  }


  if (length(alpha)!=1 || length(beta)!=1) {
    # Simulate new counts from Beta-Binomial
    # sim_y <- posterior_predict(m1.prawn.a, ndraws = 4000)  # From BRMS BB model
    S <- length(alpha_hat)
    # Initialize matrix to store simulations
    sim_y <- matrix(nrow = D, ncol = S)

    # Loop through each posterior sample and generate D Beta-Binomial samples
    for (s in 1:S) {

      if (approximate.model) {
        probs <- rbeta(D, alpha_hat[s], beta_hat[s]/group_size)  # Generate D probabilities
        sim_y[, s] <- rbinom(D, size = size, prob = probs)  # Generate D counts
      } else {
        probs <- rbeta(D, alpha_hat[s], beta_hat[s])  # Generate D probabilities
        probs_m <- 1 - (1-probs)^group_size
        sim_y[, s] <- rbinom(D, size = size, prob = probs_m)  # Generate D counts
      }
    }

    F_ymin1 <- rowMeans((sim_y) < (ObservedResponse))  # Empirical CDF
    F_y <- rowMeans((sim_y) <= ObservedResponse)

  }

  # Generate Dunn-Smyth residuals : For all domains
  epsilon <- 1e-12
  F_ymin1 <- pmin(pmax(F_ymin1, epsilon), 1 - epsilon)
  F_y <- pmin(pmax(F_y, epsilon), 1 - epsilon)
  DS <- runif(D, F_ymin1, F_y)
  normalized_DS <- qnorm(pmin(pmax(DS, 1e-10), 1 - 1e-10))
  KS <- ks.test(DS, "punif")  # Re-run KS test


  # Create a data frame for plotting
  qq_data <- data.frame(
    residuals = normalized_DS,
    group = ifelse(ObservedResponse > 0, "Non-zero", "Zero")  # Define group labels
  )

  # Perform normal Q-Q calculations
  qq <- qqnorm(qq_data$residuals, plot.it = FALSE)  # Get theoretical and sample quantiles

  # Convert into a dataframe for ggplot
  qq_df <- data.frame(
    theoretical = qq$x,
    sample = qq$y,
    group = qq_data$group
  )

  # KS test result
  ks_p_value <- round(KS$p.value, 4)

  # Plot using ggplot2
  DS_QQ_norm <- ggplot(qq_df, aes(x = theoretical, y = sample, color = group)) +
    geom_point(alpha = 0.7) +  # Q-Q points
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1.2) +  # Reference line
    labs(
      title = "Normal Q-Q Plot: Dunn-Smyth scaled residuals",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles",
      color = "Group"
    ) +
    scale_color_manual(values = c("Zero" = "grey", "Non-zero" = "blue")) +  # Define colors
    annotate("text", x = min(qq_df$theoretical), y = max(qq_df$sample),
             label = paste("KS test for Uniform: p =", ks_p_value),
             hjust = 0, vjust = 1, size = 4) +  # KS test annotation
    theme_minimal()  +
    theme(
      plot.title = element_text(size = 10),
      legend.position = "bottom",   # Move legend to bottom
      legend.title = element_text(size = 9),  # Adjust legend title size
      legend.text = element_text(size = 8)    # Adjust legend text size
    )

  # Create a data frame for plotting
  df <- data.frame(residuals = DS)
  # Plot histogram using ggplot2
  breaks <- seq(0, 1, length.out = 31)  # 30 bins over [0,1]
  # hist(df$residuals, breaks = breaks)

  DS_hist <- ggplot(df, aes(x = residuals)) +
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

  list(DS=DS,normalized_DS=normalized_DS,KS=KS,DS_hist=DS_hist,DS_QQ_norm=DS_QQ_norm)

}
DunnSmythTestBB_Custom <- function(ObservedResponse,size,group_size,alpha,beta,approximate.model=TRUE) {
  # ObservedResponse - Vector of dimension D
  # SimulatedResponse - Matrix of Dimension D X S
  # D for the length of observations
  # S for the number of simulated response
  # group_size for pool size
  # approximate.model for approximate distribution of Phi_i ~ Beta (alpha,beta/group_size)


  library(statmod)  # For qnorm and randomization
  library(VGAM)     # For pbetabinom.ab function
  library(ggplot2)

  ObservedResponse <- sort(ObservedResponse)

  D <- length(ObservedResponse) # Number of observed points
  # S <- dim(SimulatedResponse)[2]
  if (length(alpha)==1) {
    alpha_hat = alpha
  } else {alpha_hat=as.vector(alpha)}

  if (length(beta)==1) {
    beta_hat = beta
  } else {beta_hat=as.vector(beta)}

  if (length(alpha)==1 || length(beta)==1) {
    size_vector <- rep(size,D)
    alpha_hat_vector <- rep(alpha_hat,D)
    beta_hat_vector <- rep(beta_hat,D)
    valid_idx <- ObservedResponse > 0  # Identify indices where y_obs > 0
    F_ymin1 <- rep(0,length(ObservedResponse))

    if (approximate.model) {
      F_ymin1[valid_idx] <- pbetabinom.ab(ObservedResponse[valid_idx] - 1, size_vector[valid_idx], alpha_hat_vector[valid_idx], beta_hat_vector[valid_idx]/group_size)
      F_y <- pbetabinom.ab(ObservedResponse, size_vector, alpha_hat_vector, beta_hat_vector/group_size)
    } else {

      # probs <- rbeta(D,unique(alpha_hat_vector), unique(beta_hat_vector))
      # probs_m <- 1 - (1-probs)^group_size

      N_samples <- 1000
      probs <- rowMeans(matrix(rbeta(D * N_samples, unique(alpha_hat_vector), unique(beta_hat_vector)), ncol = N_samples))
      probs_m <- 1 - (1 - probs)^group_size

      F_ymin1[valid_idx] <- pbinom(ObservedResponse[valid_idx] - 1, size_vector[valid_idx], probs_m[valid_idx])
      F_y <- pbinom(ObservedResponse, size_vector, probs_m)
    }
  }


  if (length(alpha)!=1 || length(beta)!=1) {
    # Simulate new counts from Beta-Binomial
    # sim_y <- posterior_predict(m1.prawn.a, ndraws = 4000)  # From BRMS BB model
    S <- length(alpha_hat)
    # Initialize matrix to store simulations
    sim_y <- matrix(nrow = D, ncol = S)

    # Loop through each posterior sample and generate D Beta-Binomial samples
    for (s in 1:S) {

      if (approximate.model) {
        probs <- rbeta(D, alpha_hat[s], beta_hat[s]/group_size)  # Generate D probabilities
        sim_y[, s] <- rbinom(D, size = size, prob = probs)  # Generate D counts
      } else {
        probs <- rbeta(D, alpha_hat[s], beta_hat[s])  # Generate D probabilities
        probs_m <- 1 - (1-probs)^group_size
        sim_y[, s] <- rbinom(D, size = size, prob = probs_m)  # Generate D counts
      }
    }

    F_ymin1 <- rowMeans((sim_y) < (ObservedResponse))  # Empirical CDF
    F_y <- rowMeans((sim_y) <= ObservedResponse)

  }

  # Generate Dunn-Smyth residuals : For all domains
  epsilon <- 1e-12
  F_ymin1 <- pmin(pmax(F_ymin1, epsilon), 1 - epsilon)
  F_y <- pmin(pmax(F_y, epsilon), 1 - epsilon)
  DS <- runif(D, F_ymin1, F_y)
  normalized_DS <- qnorm(pmin(pmax(DS, 1e-10), 1 - 1e-10))
  KS <- ks.test(DS, "punif")  # Re-run KS test


  # Create a data frame for plotting
  qq_data <- data.frame(
    residuals = normalized_DS,
    group = ifelse(ObservedResponse > 0, "Non-zero", "Zero")  # Define group labels
  )

  # Perform normal Q-Q calculations
  qq <- qqnorm(qq_data$residuals, plot.it = FALSE)  # Get theoretical and sample quantiles

  # Convert into a dataframe for ggplot
  qq_df <- data.frame(
    theoretical = qq$x,
    sample = qq$y,
    group = qq_data$group
  )

  # KS test result
  ks_p_value <- round(KS$p.value, 4)

  # Plot using ggplot2
  DS_QQ_norm <- ggplot(qq_df, aes(x = theoretical, y = sample, color = group)) +
    geom_point(alpha = 0.7) +  # Q-Q points
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", size = 1.2) +  # Reference line
    labs(
     # title = "Normal Q-Q Plot: Dunn-Smyth scaled residuals",
      x = "Theoretical Quantiles",
      y = "Sample Quantiles",
      color = "Group"
    ) +
    scale_color_manual(values = c("Zero" = "grey", "Non-zero" = "blue")) +  # Define colors
 #   annotate("text", x = min(qq_df$theoretical), y = max(qq_df$sample),
 #            label = paste("KS test for Uniform: p =", ks_p_value),
 #            hjust = 0, vjust = 1, size = 4) +  # KS test annotation
    theme_minimal()  +
    theme(
      plot.title = element_text(size = 10),
      legend.position = "bottom",   # Move legend to bottom
      legend.title = element_blank(),  # Adjust legend title size
      legend.text = element_text(size = 8)    # Adjust legend text size
    )

  # Create a data frame for plotting
  df <- data.frame(residuals = DS)
  # Plot histogram using ggplot2
  breaks <- seq(0, 1, length.out = 31)  # 30 bins over [0,1]
  # hist(df$residuals, breaks = breaks)

  DS_hist <- ggplot(df, aes(x = residuals)) +
    geom_histogram(breaks = breaks, fill = "lightblue", color = "black", alpha = 0.7) +
    labs(
   #   title = "Distribution of Dunn-Smyth Residuals",
      x = "Residuals",
      y = "Frequency"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 10),
      axis.title = element_text(size = 9),
      axis.text = element_text(size = 8)
    )


  list(DS=DS,normalized_DS=normalized_DS,KS=KS,DS_hist=DS_hist,DS_QQ_norm=DS_QQ_norm)

}


#-------------------------------------------------------------
# Seafood Data
#-------------------------------------------------------------@

df.seafood <- read.csv("JASA Submission/deidentified_frozen_seafood.csv",header = T)
df.seafood$Yfac <- factor(df.seafood$numPos,levels=c(0:13))
x <- as.numeric(names(table(df.seafood$Yfac)))
freq <- as.numeric(table(df.seafood$Yfac))
size <- 13 # Number of bags (b)
seafood.data <- data.frame(ty=df.seafood$Yfac,n=rep(size,length(df.seafood$Yfac)))
seafood.data$ID <- c(1:dim(seafood.data)[1])
seafood.data$ty <- as.numeric(paste(seafood.data$ty))
seafood.data <- seafood.data[order(seafood.data$ty),]


# ------------------------------------------------------------------#
# Dunn-Smyth Residual Goodness-of-Fit: MH Estimator with known Imperfect Test Sensitivity (0.80)
# ------------------------------------------------------------------#

# Load MH-fitted model output (imperfect test sensitivity = 0.80)
load("JASA Submission/MH.alpha.mu.sigma.20.known.sn.80.Prawn.Rdata")

# Extract observed response (sorted)
ObservedResponse <- sort(seafood.data$ty)

# Run Dunn-Smyth residual diagnostic using posterior samples
set.seed(123456)
DSTest_MH_Imperfect <- DunnSmythTestBB_Custom(
  ObservedResponse = ObservedResponse,
  size = 13,
  group_size = 5,
  approximate.model = FALSE,
  alpha = unlist(MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$alpha_sample),
  beta = unlist(MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$beta_sample)
)

# Display diagnostic plots
DSTest_MH_Imperfect$DS_hist       # Histogram of Dunn-Smyth residuals
DSTest_MH_Imperfect$DS_QQ_norm    # Normal Q-Q plot of scaled residuals

# Print KS test result for Uniformity
DSTest_MH_Imperfect$KS

# Optional base R checks
hist(DSTest_MH_Imperfect$DS, main = "Histogram of DS Residuals")
qqnorm(DSTest_MH_Imperfect$normalized_DS)
qqline(DSTest_MH_Imperfect$normalized_DS)

# Combine plots and export
library(ggpubr)
annotate_figure(
  ggarrange(
    DSTest_MH_Imperfect$DS_hist,
    DSTest_MH_Imperfect$DS_QQ_norm,
    ncol = 2, nrow = 1
  )
)
ggsave("Dunn_Smyth_Test_MH_Estimator_imperfect_Case_sn_0.80.png", height = 4, width = 8, units = "in")



# ------------------------------------------------------------------#
# Dunn-Smyth Residual Goodness-of-Fit: MH Estimator with Unknown Test Sensitivity ~ U(0.75, 0.80)
# ------------------------------------------------------------------#

# Load MH-fitted model output with sensitivity ~ Uniform(0.75, 0.80)
load("JASA Submission/MH.alpha.mu.sigma.20.unknown.sn.80.Prawn.Rdata")

# Extract observed response (sorted)
ObservedResponse <- sort(seafood.data$ty)

# Run Dunn-Smyth residual diagnostic using posterior samples
set.seed(123456)
DSTest_MH_Imperfect_Uniform <- DunnSmythTestBB_Custom(
  ObservedResponse = ObservedResponse,
  size = 13,
  group_size = 5,
  approximate.model = FALSE,
  alpha = unlist(MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$alpha_sample),
  beta = unlist(MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$beta_sample)
)

# Display diagnostic plots
DSTest_MH_Imperfect_Uniform$DS_hist       # Histogram of Dunn-Smyth residuals
DSTest_MH_Imperfect_Uniform$DS_QQ_norm    # Normal Q-Q plot of scaled residuals

# Print KS test result for Uniformity
DSTest_MH_Imperfect_Uniform$KS

# Optional base R diagnostic plots
hist(DSTest_MH_Imperfect_Uniform$DS, main = "Histogram of DS Residuals")
qqnorm(DSTest_MH_Imperfect_Uniform$normalized_DS)
qqline(DSTest_MH_Imperfect_Uniform$normalized_DS)

# Combine ggplots and export to file
annotate_figure(
  ggarrange(
    DSTest_MH_Imperfect_Uniform$DS_hist,
    DSTest_MH_Imperfect_Uniform$DS_QQ_norm,
    ncol = 2, nrow = 1
  )
)
ggsave("Dunn_Smyth_Test_MH_Estimator_imperfect_Case_sn_Unif_0.80.png", height = 4, width = 8, units = "in")



