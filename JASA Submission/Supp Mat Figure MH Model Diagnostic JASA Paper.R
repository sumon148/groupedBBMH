# ============================================================================ #
# Title    :  Diagnostics of MCMC outputs under Metropolis-Hastings Algorithm
# Project  : JASA Submission
# ============================================================================ #

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

# --------------------------------------------------------------------------- #
# Load MCMC Results for Known and Unknown Sensitivity Models
# --------------------------------------------------------------------------- #
load("JASA Submission/MH.alpha.mu.sigma.20.known.sn.80.Prawn.Rdata")
load("JASA Submission/MH.alpha.mu.sigma.20.unknown.sn.80.Prawn.Rdata")

# --------------------------------------------------------------------------- #
# Extract and Prepare Parameters for Diagnostic Plots
# --------------------------------------------------------------------------- #
selected_params_known_imPerfect_sensitivity <- lapply(MH.alpha.mu.sigma.20.known.sn.80.Prawn$parameters, function(mat) mat[, 1:2])
var.names <- c("log(alpha)","logit(mu)")
prawn_mcmc_known_imperfect_sensitivity <- create_mcmc_BP(selected_params_known_imPerfect_sensitivity,var.names)

selected_params_unknown_imPerfect_sensitivity <- lapply(MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$parameters, function(mat) mat[, 1:3])
var.names <- c("log(alpha)","logit(mu)","logit(Sensitivity)")
prawn_mcmc_unknown_imperfect_sensitivity <- create_mcmc_BP(selected_params_unknown_imPerfect_sensitivity,var.names)


# --------------------------------------------------------------------------- #
# Set color scheme for plots
# --------------------------------------------------------------------------- #
color_scheme_set("viridis")

# --------------------------------------------------------------------------- #
# Traceplots for Known and Unknown Sensitivity Models
# --------------------------------------------------------------------------- #
p1_Prawn <- mcmc_trace(prawn_mcmc_known_imperfect_sensitivity) # Trace plot for known sensitivity
p2_Prawn <- mcmc_trace(prawn_mcmc_unknown_imperfect_sensitivity) # Trace plot for unknown sensitivity

# --------------------------------------------------------------------------- #
# Density Overlay Plots for Known and Unknown Sensitivity Models
# --------------------------------------------------------------------------- #
p1_Prawn_dens_overlay <- mcmc_dens_overlay(prawn_mcmc_known_imperfect_sensitivity)  # Density overlay for known sensitivity
p2_Prawn_dens_overlay <- mcmc_dens_overlay(prawn_mcmc_unknown_imperfect_sensitivity) # Density overlay for unknown sensitivity

# --------------------------------------------------------------------------- #
# Arrange plots vertically and save to file
# --------------------------------------------------------------------------- #

# Combine and save known sensitivity plots
ggarrange(p1_Prawn, p1_Prawn_dens_overlay, ncol = 1, common.legend = TRUE)
# ggsave("Trace_Plot_Prawn_known_sensitivity.png", width = 6, height = 8)

# Combine and save unknown sensitivity plots
ggarrange(p2_Prawn, p2_Prawn_dens_overlay, ncol = 1, common.legend = TRUE)
#ggsave("Trace_Plot_Prawn_unknown_sensitivity.png", width = 8, height = 8)


# --------------------------------------------------------------------------- #
# Leakage Estimate for Prawn Data: Known and Unknown Imperfect Sensitivity
# --------------------------------------------------------------------------- #

# Estimate leakage samples for known imperfect sensitivity (sensitivity fixed at 0.80)
prawn.Leakage.post.samples.known.imprf.sn <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                                    alpha=MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$alpha_sample,
                                                    beta=MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$beta_sample,
                                                    sensitivity = 0.80,
                                                    cutoff = 0)
# Estimate leakage samples for unknown imperfect sensitivity (sensitivity follows U(0.75,0.85))
prawn.Leakage.post.samples.unknown.imprf.sn <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                                          alpha=MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$alpha_sample,
                                                          beta=MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$beta_sample,
                                                          sensitivity=MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$sensitivity_sample,
                                                          cutoff = 0)
# --------------------------------------------------------------------------- #
# Summary statistics for MCMC samples - Known Imperfect Sensitivity
# --------------------------------------------------------------------------- #
prawn.post.summ.known.imprf.sn <- summary_mcmc(
  MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$alpha_sample,
  MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$beta_sample,
  MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$mu_sample,
  prawn.Leakage.post.samples.known.imprf.sn$E.L,
  prawn.Leakage.post.samples.known.imprf.sn$P.L,
  varnames = c("alpha","beta","mu","E.L","P.L"))
Rhat.post.summ.known.imprf <- prawn.post.summ.known.imprf.sn$summary_mcmc[c("alpha","beta","mu","E.L","P.L"),c("R_hat")]

# --------------------------------------------------------------------------- #
# Summary statistics for MCMC samples - Unknown Imperfect Sensitivity
# --------------------------------------------------------------------------- #

prawn.post.summ.unknown.imprf.sn <- summary_mcmc(
  MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$alpha_sample,
  MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$beta_sample,
  MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$mu_sample,
  MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$sensitivity_sample,
  prawn.Leakage.post.samples.unknown.imprf.sn$E.L,
  prawn.Leakage.post.samples.unknown.imprf.sn$P.L,
  varnames = c("alpha","beta","mu","sensitivity","E.L","P.L")) #

Rhat.post.summ.unknown.imprf <- prawn.post.summ.unknown.imprf.sn$summary_mcmc[c("alpha","beta","mu","E.L","P.L"),c("R_hat")]

# --------------------------------------------------------------------------- #
# Combine summary statistics for known and unknown imperfect sensitivity
# --------------------------------------------------------------------------- #
summary_stats_mcmc <- rbind(
  prawn.post.summ.known.imprf.sn$summary_mcmc,
  prawn.post.summ.unknown.imprf.sn$summary_mcmc
)

library(xtable)
xtable(summary_stats_mcmc,digits = 3)


