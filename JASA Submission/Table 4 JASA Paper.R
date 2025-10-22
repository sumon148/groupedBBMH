
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

# ----------------------------------------------------------
# Table 4: Estimation of Beta and Leakage Parameters
#          by Minimum Positive Propensity
# ----------------------------------------------------------

# Read the cleaned and de-identified dataset
df_prawn <- read.csv("JASA Submission/deidentified_frozen_seafood.csv")

# Define outcome as a factor based on number of positive samples (0 to 13)
df_prawn$Yfac <- factor(df_prawn$numPos, levels = 0:13)

# Extract unique outcome levels and their frequencies
outcome_levels <- as.numeric(names(table(df_prawn$Yfac)))
frequencies <- as.numeric(table(df_prawn$Yfac))

# Define trial size (number of samples per unit)
trial_size <- 13  # Number of bags per sample

# Construct binomial response data frame
prawn_data <- data.frame(
  ty = as.numeric(as.character(df_prawn$Yfac)),
  n = rep(trial_size, nrow(df_prawn)),
  ID = 1:nrow(df_prawn)
)

# Sort data by outcome value
prawn_data <- prawn_data[order(prawn_data$ty), ]

# Create summary table by unique outcome value
summary_prawn_data <- as.data.frame(table(prawn_data$ty))
colnames(summary_prawn_data) <- c("ty", "freq")
summary_prawn_data$ty <- as.numeric(as.character(summary_prawn_data$ty))

# ----------------------------------------------------------
# Define input variables for the MH (Metropolis-Hastings) algorithm
# ----------------------------------------------------------

b <- trial_size             # Number of trials per observation
Nbar <- 5                   # Prior or hyperparameter (as required)
ty <- prawn_data$ty         # Observed positive counts
D <- length(ty)             # Total number of observations

ty2 <- outcome_levels       # Unique positive counts
wt <- frequencies           # Frequency of each unique outcome
d <- length(ty2)            # Number of unique outcome groups

# -------------------------------------------------------------------------
# Apply MH Algorithm to Fit Generalized Beta-Binomial Models
# Settings:
#   - Delta = 0.80 (fixed sensitivity)
#   - Lambda = 1.0
#   - G = 5000 (number of quadrature points)
#   - R = 50,000 (MCMC iterations per chain)
#   - Burn-in = 50%
#   - Thinning = every 5th draw
#   - 4 chains with fixed seeds
#   - cutoff (k) values: 0%, 0.25%, 0.5%, 1%, 2%
# -------------------------------------------------------------------------

# k=0.0000

G=5000;R=50e3;
par <- c(log(0.005), logit(7e-04), logit(0.70)) ; initial <- par; k <- 0.0000
MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0 <- MH.Sampler.Chain.TBB(par=par,cutoff=k,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                                      sigma=0.2,num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=TRUE,specificity=FALSE,
                                                                      sensitivity.range=c(0.80,0.80),specificity.range=NULL)
save(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0,file="JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0.Rdata")


# k=0.0025

G=5000;R=50e3;
par <- c(log(0.005), logit(7e-04), logit(0.70)) ; initial <- par; k <- 0.0025
MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a <- MH.Sampler.Chain.TBB(par=par,cutoff=k,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                                      sigma=0.2,num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=TRUE,specificity=FALSE,
                                                                      sensitivity.range=c(0.80,0.80),specificity.range=NULL)
save(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a,file="JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a.Rdata")


# k=0.005
G=5000;R=50e3;
par <- c(log(0.005), logit(7e-04), logit(0.70)) ; initial <- par; k <- 0.005
MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b <- MH.Sampler.Chain.TBB(par=par,cutoff=k,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                                      sigma=0.2,num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=TRUE,specificity=FALSE,
                                                                      sensitivity.range=c(0.80,0.80),specificity.range=NULL)
save(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b,file="JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b.Rdata")

# k=0.01
G=5000;R=50e3;
par <- c(log(0.005), logit(7e-04), logit(0.70)) ; initial <- par; k <- 0.010
MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c <- MH.Sampler.Chain.TBB(par=par,cutoff=k,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                                      sigma=0.2,num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=TRUE,specificity=FALSE,
                                                                      sensitivity.range=c(0.80,0.80),specificity.range=NULL)
save(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c,file="JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c.Rdata")

# k=0.02
G=5000;R=50e3;
par <- c(log(0.005), logit(7e-04), logit(0.70)) ; initial <- par; k <- 0.020
MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d <- MH.Sampler.Chain.TBB(par=par,cutoff=k,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                                      sigma=0.2,num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=TRUE,specificity=FALSE,
                                                                      sensitivity.range=c(0.80,0.80),specificity.range=NULL)
save(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d,file="JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d.Rdata")



# ---------------------------------------------------------#
# Load MCMC outputs of the fitted Generalized BB Models----
# Using Delta=0.80 as imperfect situation
# k=0.0025, 0.0050, 0.01, 0.02
# -------------------------------------------------#

load("JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0.Rdata") # k=0
load("JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a.Rdata") # k=0.0025
load("JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b.Rdata") # k=0.0050
load("JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c.Rdata") # k=0.0100
load("JASA Submission/MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d.Rdata") # k=0.0200
MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0 <- MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0
MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.a <- MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a
MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.b <- MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b
MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c <- MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c
MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d <- MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d


# ---------------------------------------------------------#
# Extract and summarize Beta Parameters -------
# ----------------------------------------------------------#

summ.alpha.beta.TBB.00 <- summary_mcmc(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$alpha_sample,
                                       MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$beta_sample,
                                       MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$mu_sample,
                                       varnames = c("alpha","beta","mu"))
summ.alpha.beta.TBB.00$summary_mcmc

summ.alpha.beta.TBB.0025 <- summary_mcmc(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.a$target.parameters$alpha_sample,
                                         MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.a$target.parameters$beta_sample,
                                         MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.a$target.parameters$mu_sample,
                                         varnames = c("alpha","beta","mu"))
summ.alpha.beta.TBB.0025$summary_mcmc

summ.alpha.beta.TBB.005 <- summary_mcmc(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.b$target.parameters$alpha_sample,
                                        MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.b$target.parameters$beta_sample,
                                        MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.b$target.parameters$mu_sample,
                                        varnames = c("alpha","beta","mu"))
summ.alpha.beta.TBB.005$summary_mcmc

summ.alpha.beta.TBB.01 <- summary_mcmc(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$alpha_sample,
                                       MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$beta_sample,
                                       MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$mu_sample,
                                       varnames = c("alpha","beta","mu"))
summ.alpha.beta.TBB.01$summary_mcmc

summ.alpha.beta.TBB.02 <- summary_mcmc(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$alpha_sample,
                                       MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$beta_sample,
                                       MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$mu_sample,
                                       varnames = c("alpha","beta","mu"))
summ.alpha.beta.TBB.02$summary_mcmc


# ---------------------------------------------------------#
# Leakage estimates with Sensitivity=0.8 -------
# ----------------------------------------------------------#

E.P.trBB.cutoff.00.imperfect.80 <- estimate_leakage(alpha = unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$alpha_sample),
                                                    beta = unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$beta_sample),
                                                    cutoff = 0,
                                                    B = 8000,M = 40,b = 13,m = 5,
                                                    sensitivity = 0.80)
E.P.trBB.cutoff.01.imperfect.80 <- estimate_leakage(alpha = unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$alpha_sample),
                                                    beta = unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$alpha_sample),
                                                    cutoff = 0.01,
                                                    B = 8000,M = 40,b = 13,m = 5,
                                                    sensitivity = 0.80)
E.P.trBB.cutoff.02.imperfect.80 <- estimate_leakage(alpha = unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$alpha_sample),
                                                    beta = unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$alpha_sample),
                                                    cutoff = 0.02,
                                                    B = 8000,M = 40,b = 13,m = 5,
                                                    sensitivity = 0.80)


# ----------------------------------------------------------------------
# Prepare Table 4: Posterior Summaries for Parameters at Different Cutoffs
# ----------------------------------------------------------------------

summarise_Mean_CI <- function(x) {
  c(mean = mean(x, na.rm = TRUE),
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

Summ.alpha.cutoff.00.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$alpha_sample))*10^3
Summ.alpha.cutoff.01.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$alpha_sample))*10^3
Summ.alpha.cutoff.02.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$alpha_sample))*10^3

Summ.beta.cutoff.00.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$beta_sample))
Summ.beta.cutoff.01.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$beta_sample))
Summ.beta.cutoff.02.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$beta_sample))

Summ.mu.cutoff.00.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$mu_sample))*10^3
Summ.mu.cutoff.01.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$mu_sample))*10^3
Summ.mu.cutoff.02.imperfect.80 <- summarise_Mean_CI(unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$mu_sample))*10^3

Summ.E.L.cutoff.00.imperfect.80 <- summarise_Mean_CI(E.P.trBB.cutoff.00.imperfect.80$E.L)
Summ.E.L.cutoff.01.imperfect.80 <- summarise_Mean_CI(E.P.trBB.cutoff.01.imperfect.80$E.L)
Summ.E.L.cutoff.02.imperfect.80 <- summarise_Mean_CI(E.P.trBB.cutoff.02.imperfect.80$E.L)
Summ.P.L.cutoff.00.imperfect.80 <- summarise_Mean_CI(E.P.trBB.cutoff.00.imperfect.80$P.L)*10^2
Summ.P.L.cutoff.01.imperfect.80 <- summarise_Mean_CI(E.P.trBB.cutoff.01.imperfect.80$P.L)*10^2
Summ.P.L.cutoff.02.imperfect.80 <- summarise_Mean_CI(E.P.trBB.cutoff.02.imperfect.80$P.L)*10^2

Table_4 <- as.data.frame(
  rbind(
c(Summ.alpha.cutoff.00.imperfect.80,Summ.alpha.cutoff.01.imperfect.80,Summ.alpha.cutoff.02.imperfect.80),
c(Summ.beta.cutoff.00.imperfect.80,Summ.beta.cutoff.01.imperfect.80,Summ.beta.cutoff.02.imperfect.80),
c(Summ.mu.cutoff.00.imperfect.80,Summ.mu.cutoff.01.imperfect.80,Summ.mu.cutoff.02.imperfect.80),
c(Summ.E.L.cutoff.00.imperfect.80,Summ.E.L.cutoff.01.imperfect.80,Summ.E.L.cutoff.02.imperfect.80),
c(Summ.P.L.cutoff.00.imperfect.80,Summ.P.L.cutoff.01.imperfect.80,Summ.P.L.cutoff.02.imperfect.80)
)
)
rownames(Table_4) <- c("alpha*10^3","beta","mu*10^3","E.L","P.L*10^2")
colnames(Table_4) <- c("k0:Mean","k0:Q_2.5%","k0:Q_97.5%","k01:Mean","k01:Q_2.5%","k01:Q_97.5%","k02:Mean","k02:Q_2.5%","k02:Q_97.5%")

library(xtable)
xtable(Table_4,digits=2)


# --------------------------------------------------------------------------
# Supplementary Table: Leakage Estimates under Varying Sensitivity
# Assumption: Perfect Specificity (i.e., Specificity = 1.00)
# Goal: Evaluate expected leakage across different sensitivity scenarios
# --------------------------------------------------------------------------

alpha_tbb_00 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$alpha_sample)
beta_tbb_00 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$beta_sample)
mu_tbb_00 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.0$target.parameters$mu_sample)
alpha_tbb_0025 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.a$target.parameters$alpha_sample)
beta_tbb_0025 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.a$target.parameters$beta_sample)
mu_tbb_0025 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.a$target.parameters$mu_sample)
alpha_tbb_005 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.b$target.parameters$alpha_sample)
beta_tbb_005 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.b$target.parameters$beta_sample)
mu_tbb_005 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.b$target.parameters$mu_sample)
alpha_tbb_01 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$alpha_sample)
beta_tbb_01 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$beta_sample)
mu_tbb_01 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.c$target.parameters$mu_sample)
alpha_tbb_02 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$alpha_sample)
beta_tbb_02 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$beta_sample)
mu_tbb_02 <- unlist(MH.alpha.mu.sigma.20.Imperfect.Prawn.TBB.d$target.parameters$mu_sample)

E.P.trBB.cutoff.00.perfect <- estimate_leakage(alpha = alpha_tbb_00,beta = beta_tbb_00,cutoff = 0,B = 8000,M = 40,b = 13,m = 5,sensitivity = 1)
E.P.trBB.cutoff.00.imperfect.95 <- estimate_leakage(alpha = alpha_tbb_00,beta = beta_tbb_00,cutoff = 0,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.95)
E.P.trBB.cutoff.00.imperfect.90 <- estimate_leakage(alpha = alpha_tbb_00,beta = beta_tbb_00,cutoff = 0,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.90)
E.P.trBB.cutoff.00.imperfect.80 <- estimate_leakage(alpha = alpha_tbb_00,beta = beta_tbb_00,cutoff = 0,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.80)
E.P.trBB.cutoff.00.imperfect.70 <- estimate_leakage(alpha = alpha_tbb_00,beta = beta_tbb_00,cutoff = 0,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.70)

E.P.trBB.cutoff.0025.perfect <- estimate_leakage(alpha = alpha_tbb_0025,beta = beta_tbb_0025,cutoff = 0.0025,B = 8000,M = 40,b = 13,m = 5,sensitivity = 1)
E.P.trBB.cutoff.0025.imperfect.95 <- estimate_leakage(alpha = alpha_tbb_0025,beta = beta_tbb_0025,cutoff = 0.0025,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.95)
E.P.trBB.cutoff.0025.imperfect.90 <- estimate_leakage(alpha = alpha_tbb_0025,beta = beta_tbb_0025,cutoff = 0.0025,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.90)
E.P.trBB.cutoff.0025.imperfect.80 <- estimate_leakage(alpha = alpha_tbb_0025,beta = beta_tbb_0025,cutoff = 0.0025,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.80)
E.P.trBB.cutoff.0025.imperfect.70 <- estimate_leakage(alpha = alpha_tbb_0025,beta = beta_tbb_0025,cutoff = 0.0025,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.70)

E.P.trBB.cutoff.005.perfect <- estimate_leakage(alpha = alpha_tbb_005,beta = beta_tbb_005,cutoff = 0.005,B = 8000,M = 40,b = 13,m = 5,sensitivity = 1)
E.P.trBB.cutoff.005.imperfect.95 <- estimate_leakage(alpha = alpha_tbb_005,beta = beta_tbb_005,cutoff = 0.005,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.95)
E.P.trBB.cutoff.005.imperfect.90 <- estimate_leakage(alpha = alpha_tbb_005,beta = beta_tbb_005,cutoff = 0.005,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.90)
E.P.trBB.cutoff.005.imperfect.80 <- estimate_leakage(alpha = alpha_tbb_005,beta = beta_tbb_005,cutoff = 0.005,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.80)
E.P.trBB.cutoff.005.imperfect.70 <- estimate_leakage(alpha = alpha_tbb_005,beta = beta_tbb_005,cutoff = 0.005,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.70)


E.P.trBB.cutoff.01.perfect <- estimate_leakage(alpha = alpha_tbb_01,beta = beta_tbb_01,cutoff = 0.01,B = 8000,M = 40,b = 13,m = 5,sensitivity = 1)
E.P.trBB.cutoff.01.imperfect.95 <- estimate_leakage(alpha = alpha_tbb_01,beta = beta_tbb_01,cutoff = 0.01,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.95)
E.P.trBB.cutoff.01.imperfect.90 <- estimate_leakage(alpha = alpha_tbb_01,beta = beta_tbb_01,cutoff = 0.01,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.90)
E.P.trBB.cutoff.01.imperfect.80 <- estimate_leakage(alpha = alpha_tbb_01,beta = beta_tbb_01,cutoff = 0.01,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.80)
E.P.trBB.cutoff.01.imperfect.70 <- estimate_leakage(alpha = alpha_tbb_01,beta = beta_tbb_01,cutoff = 0.01,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.70)


E.P.trBB.cutoff.02.perfect <- estimate_leakage(alpha = alpha_tbb_02,beta = beta_tbb_02,cutoff = 0.02,B = 8000,M = 40,b = 13,m = 5,sensitivity = 1)
E.P.trBB.cutoff.02.imperfect.95 <- estimate_leakage(alpha = alpha_tbb_02,beta = beta_tbb_02,cutoff = 0.02,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.95)
E.P.trBB.cutoff.02.imperfect.90 <- estimate_leakage(alpha = alpha_tbb_02,beta = beta_tbb_02,cutoff = 0.02,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.90)
E.P.trBB.cutoff.02.imperfect.80 <- estimate_leakage(alpha = alpha_tbb_02,beta = beta_tbb_02,cutoff = 0.02,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.80)
E.P.trBB.cutoff.02.imperfect.70 <- estimate_leakage(alpha = alpha_tbb_02,beta = beta_tbb_02,cutoff = 0.02,B = 8000,M = 40,b = 13,m = 5,sensitivity = 0.70)


# --------------------------------------------------------------------------- #
# Step 1: Prepare Leakage Data by Sensitivity and Cutoff Thresholds
# Description: Converts output lists to data frames for comparative analysis
# across different combinations of test sensitivity and cutoff values.
# --------------------------------------------------------------------------- #

convert_list_to_df <- function(obj, cutoff, sensitivity) {
  n <- length(obj$alpha)  # assume all vector elements are of same length

  data.frame(
    alpha      = obj$alpha,
    beta       = obj$beta,
    mu         = obj$mu,
    E.L        = obj$E.L,
    P.L        = obj$P.L,
    B          = rep(obj$B, n),
    b          = rep(obj$b, n),
    M          = rep(obj$M, n),
    m          = rep(obj$m, n),
    cutoff     = rep(cutoff, n),
    sensitivity = rep(sensitivity, n)
  )
}

# --------------------------------------------------------------------------- #
# Step 2: Generate Data Frames for Each Cutoff × Sensitivity Combination
# --------------------------------------------------------------------------- #

df_00_1.0 <- convert_list_to_df(E.P.trBB.cutoff.00.perfect, cutoff=0.00, sensitivity=1.0)
df_00_0.95 <- convert_list_to_df(E.P.trBB.cutoff.00.imperfect.95,0.00, sensitivity = 0.95)
df_00_0.90 <- convert_list_to_df(E.P.trBB.cutoff.00.imperfect.90, cutoff = 0.00, sensitivity = 0.90)
df_00_0.80 <- convert_list_to_df(E.P.trBB.cutoff.00.imperfect.80, cutoff = 0.00, sensitivity = 0.80)
df_00_0.70 <- convert_list_to_df(E.P.trBB.cutoff.00.imperfect.70, cutoff = 0.00, sensitivity = 0.70)


df_0025_1.0 <- convert_list_to_df(E.P.trBB.cutoff.0025.perfect, cutoff = 0.0025, sensitivity = 1.0)
df_0025_0.95 <- convert_list_to_df(E.P.trBB.cutoff.0025.imperfect.95, cutoff = 0.0025, sensitivity = 0.95)
df_0025_0.90 <- convert_list_to_df(E.P.trBB.cutoff.0025.imperfect.90, cutoff = 0.0025, sensitivity = 0.90)
df_0025_0.80 <- convert_list_to_df(E.P.trBB.cutoff.0025.imperfect.80, cutoff = 0.0025, sensitivity = 0.80)
df_0025_0.70 <- convert_list_to_df(E.P.trBB.cutoff.0025.imperfect.70, cutoff = 0.0025, sensitivity = 0.70)


df_005_1.0 <- convert_list_to_df(E.P.trBB.cutoff.005.perfect,cutoff = 0.005, sensitivity = 1.0)
df_005_0.95 <- convert_list_to_df(E.P.trBB.cutoff.005.imperfect.95, cutoff = 0.005, sensitivity = 0.95)
df_005_0.90 <- convert_list_to_df(E.P.trBB.cutoff.005.imperfect.90, cutoff = 0.005, sensitivity = 0.90)
df_005_0.80 <- convert_list_to_df(E.P.trBB.cutoff.005.imperfect.80,cutoff = 0.005, sensitivity = 0.80)
df_005_0.70 <- convert_list_to_df(E.P.trBB.cutoff.005.imperfect.70,cutoff = 0.005, sensitivity = 0.70)


df_01_1.0 <- convert_list_to_df(E.P.trBB.cutoff.01.perfect, cutoff = 0.01, sensitivity = 1.0)
df_01_0.95 <- convert_list_to_df(E.P.trBB.cutoff.01.imperfect.95, cutoff = 0.01, sensitivity = 0.95)
df_01_0.90 <- convert_list_to_df(E.P.trBB.cutoff.01.imperfect.90, cutoff = 0.01, sensitivity = 0.90)
df_01_0.80 <- convert_list_to_df(E.P.trBB.cutoff.01.imperfect.80, cutoff = 0.01, sensitivity = 0.80)
df_01_0.70 <- convert_list_to_df(E.P.trBB.cutoff.01.imperfect.70, cutoff = 0.01, sensitivity = 0.70)


df_02_1.0 <- convert_list_to_df(E.P.trBB.cutoff.02.perfect, cutoff = 0.02, sensitivity = 1.0)
df_02_0.95 <- convert_list_to_df(E.P.trBB.cutoff.02.imperfect.95, cutoff = 0.02, sensitivity = 0.95)
df_02_0.90 <- convert_list_to_df(E.P.trBB.cutoff.02.imperfect.90, cutoff = 0.02, sensitivity = 0.90)
df_02_0.80 <- convert_list_to_df(E.P.trBB.cutoff.02.imperfect.80, cutoff = 0.02, sensitivity = 0.80)
df_02_0.70 <- convert_list_to_df(E.P.trBB.cutoff.02.imperfect.70, cutoff = 0.02, sensitivity = 0.70)


# --------------------------------------------------------------------------- #
# Step 3: Combine All Data Frames into a Single Long Format Dataset
# --------------------------------------------------------------------------- #

combined_results_Leakage_Estimate_Cutoff_Sensitivity <- rbind(
  df_00_1.0,
  df_00_0.95,
  df_00_0.90,
  df_00_0.80,
  df_00_0.70,

  df_0025_1.0,
  df_0025_0.95,
  df_0025_0.90,
  df_0025_0.80,
  df_0025_0.70,

  df_005_1.0,
  df_005_0.95,
  df_005_0.90,
  df_005_0.80,
  df_005_0.70,

  df_01_1.0,
  df_01_0.95,
  df_01_0.90,
  df_01_0.80,
  df_01_0.70,

  df_02_1.0,
  df_02_0.95,
  df_02_0.90,
  df_02_0.80,
  df_02_0.70
)


# --------------------------------------------------------------------------- #
# Step 4: Summarize Leakage Estimates by Cutoff × Sensitivity
# --------------------------------------------------------------------------- #

E.P.Li.Estimate.Cutoff.Sensitivity <- combined_results_Leakage_Estimate_Cutoff_Sensitivity %>%
  group_by(cutoff, sensitivity) %>%
  summarise(
    mean_E_Li = mean(E.L),
    lower_E_Li = quantile(E.L, probs = 0.025),
    upper_E_Li = quantile(E.L, probs = 0.975),

    mean_P_Li = mean(P.L),
    lower_P_Li = quantile(P.L, probs = 0.025),
    upper_P_Li = quantile(P.L, probs = 0.975),

    .groups = "drop"
  )

E.P.Li.Estimate.Cutoff.Sensitivity <- E.P.Li.Estimate.Cutoff.Sensitivity[order(-E.P.Li.Estimate.Cutoff.Sensitivity$sensitivity,
                                                                               E.P.Li.Estimate.Cutoff.Sensitivity$cutoff),]

colnames(E.P.Li.Estimate.Cutoff.Sensitivity) <- c("cutoff","sensitivity","E.Li","E.Li.ll",
                                                  "E.Li.ul","P.Li","P.Li.ll","P.Li.ul")

# write.csv(E.P.Li.Community.Risk,file="E.P.Li.Community.Risk.csv")
# write.csv(E.P.Li.Estimate.Cutoff.Sensitivity,file="E.P.Li.Estimate.Cutoff.Sensitivity.csv")

