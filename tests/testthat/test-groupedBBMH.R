library(testthat)
library(groupedBBMH)  #

# -----------------------------------------------------------------------------
# Load and preprocess test data
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# Shared test data setup
# -----------------------------------------------------------------------------
library(testthat)
library(groupedBBMH)

# Load the frozen seafood dataset shipped with the package
seafood_csv <- system.file("extdata", "deidentified_frozen_seafood.csv", package = "groupedBBMH")
df.seafood <- read.csv(seafood_csv)

# Convert the number of positive counts to a factor
df.seafood$Yfac <- factor(df.seafood$numPos, levels = 0:13)

# Counts of each unique value
x <- as.numeric(names(table(df.seafood$Yfac)))     # unique ty values
freq <- as.numeric(table(df.seafood$Yfac))        # frequency of each ty

# Prepare seafood.data for MH sampler tests
size <- 13  # number of bags
seafood.data <- data.frame(
  ty = as.numeric(as.character(df.seafood$Yfac)),
  n = rep(size, nrow(df.seafood)),
  ID = seq_len(nrow(df.seafood))
)
seafood.data <- seafood.data[order(seafood.data$ty), ]

# Summarized frequency table for tests
summ.seafood.data <- data.frame(
  ty = x,
  freq = freq
)

test_that("Data preprocessing works correctly", {
  expect_equal(length(seafood.data$ty), nrow(df.seafood))
  expect_equal(sum(summ.seafood.data$freq), nrow(seafood.data))
})


# -----------------------------------------------------------------------------
# MCMC / MH algorithm tests
# -----------------------------------------------------------------------------
test_that("MH_Sampler_BB and summary_mcmc run without errors", {

  b <- 13
  Nbar <- 5
  ty <- seafood.data$ty
  ty2 <- x
  wt <- freq

  G <- 500; R <- 1e3
  par <- c(log(0.005), logit(7e-04))

  MH.alpha.mu.TBB.0 <- MH_Sampler_BB(par = par, cutoff = 0,
                                     initial = par,
                                     alpha = TRUE, beta = FALSE, mu = TRUE, rho = FALSE,
                                     G = G, R = R, sigma = 0.2,
                                     num_chains = 4, burnin = R * 0.5, thin = 5,
                                     seed = c(123, 456, 789, 135),
                                     trail = b, ty = ty2, wt = wt, group.size = Nbar,
                                     sensitivity = FALSE, specificity = FALSE,
                                     sensitivity.range = NULL, specificity.range = NULL)

  summ.TBB.0 <- summary_mcmc(MH.alpha.mu.TBB.0$target.parameters$alpha_sample,
                             MH.alpha.mu.TBB.0$target.parameters$beta_sample,
                             MH.alpha.mu.TBB.0$target.parameters$mu_sample,
                             varnames = c("alpha", "beta", "mu"))

  expect_true(!is.null(summ.TBB.0$summary_mcmc))
})

# -----------------------------------------------------------------------------
# MLE / BB fitting tests
# -----------------------------------------------------------------------------
test_that("GroupedBB and truncated BB fitting works", {

  freq <- c(2815,9,10,6,1,3,2,0,1,2,1,0,0,0)
  ty <- 0:13

  bb.m1 <- fit_GroupedBB(ty = ty, freq = freq, b = 13, m = 5, sensitivity = 1, specificity = 1)
  expect_true(!is.null(bb.m1))

  bb.m2 <- fit_GroupedBB(ty = ty, freq = freq, b = 13, m = 5, M = 40, B = 8000,
                         sensitivity = 1, specificity = 1, leakage = TRUE)
  expect_true(!is.null(bb.m2))

  trbb.m1 <- fit_trGroupedBB(ty = ty, freq = freq, b = 13, m = 5, cutoff = 0)
  expect_true(!is.null(trbb.m1))

  trbb.m2 <- fit_trGroupedBB(ty, freq, b = 13, m = 5, cutoff = 0.01, sensitivity = 0.8,
                             specificity = 1, R = 10000)
  expect_true(!is.null(trbb.m2))
})

# -----------------------------------------------------------------------------
# Example MLE optimization tests
# -----------------------------------------------------------------------------
test_that("MLE optimization runs without error", {

  MLE.BB.m5 <- optim(c(log(0.005), log(3.0), log(0.5)), loglik_group_bb,
                     ty = ty, freq = freq, b = 13, theta = Inf, m = 5, M = 40,
                     sensitivity = 1, deviance = FALSE,
                     control = list(factr = 1e-12, fnscale = -1, trace = FALSE),
                     R = 1e4, hessian = FALSE, method = "L-BFGS-B",
                     lower = c(log(0.0001), log(0.0001), log(0.99)),
                     upper = c(Inf, Inf, 0))

  expect_true(!is.null(MLE.BB.m5$par))

  beta.est.cutoff.01 <- optim(c(0), loglik_group_trbb,
                              ty = ty, freq = freq, b = 13, m = 5, theta = Inf, R = 1e4,
                              alpha = 0.005, cutoff = 0.01, sensitivity = 1, specificity = 1,
                              deviance = FALSE, method = "Brent",
                              control = list(reltol = 1e-12, fnscale = -1),
                              hessian = FALSE, lower = log(0.05), upper = log(100))

  expect_true(!is.null(beta.est.cutoff.01$par))
})
