# --------------------------------------------------------------------------- #
# Load Required Libraries
# --------------------------------------------------------------------------- #
library(groupedBBMH)  # Custom package for grouped beta-binomial MH modeling
library(ggplot2)
library(ggpubr)
library(xtable)

#-------------------------------------------------------------@
# Figure 5: JASA PAper
# Supplementary Figure S5
# Supplementary Table S4
# Model-based Simulation using MH-based TrBB model under Imperfect Testing
# Community Risk Assessment
# Using Frozen-food imported data
# Exact BB Model: Bayesian method MH Algorithm
# Target Parameters: TP, TN, FN
# FN: Leakage
# P(Li>0): FN Batches
# E(Li): Number of contaminated items in FN Batches
# Delta=0.8, k=0.0, 0.0025, 0.005, 0.01, 0.02, m=5, b=13
#-------------------------------------------------------------@


# Simulations: Do Not Run-----------------------------------------------------------------------------
#
#
# Description:
#   This simulation evaluates the TrBB (Truncated Beta-Binomial) model using
#   the Metropolis-Hastings (MH) algorithm. The testing procedure is assumed
#   to be imperfect with parameters:
#     - Delta = 0.80
#     - Lambda = 1
#     - k = {0, 0.0025, 0.005, 0.01, 0.02}
#
# Posterior sample of alpha and beta values are available from:
#     The results are stored in the following RData file:
#     "MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0.Rdata"
#     "MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a.Rdata"
#     "MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b.Rdata"
#     "MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c.Rdata"
#     "MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d.Rdata"
# -----------------------------------------------------------------------------#

load("submission\\MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0.Rdata")

alpha_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0$target.parameters$alpha_sample),5000)
beta_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0$target.parameters$beta_sample),5000)

# alpha_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0$target.parameters$alpha_sample),100)
# beta_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0$target.parameters$beta_sample),100)

cut_off <- 0.00
hist(alpha_TrBB,breaks = 100)
hist(beta_TrBB,breaks = 100)

summary(alpha_TrBB)
summary(beta_TrBB)

Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 10000,B = 8000,M = 40,m = 5,b = 13,sensitivity_values = c(0.7,0.8,0.9,0.95,1.0),NSim = length(alpha_TrBB),seed = 2015)
# Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 100,B = 8000,M = 40,m = 5,b = 13,sensitivity_values = c(0.8,1.0),NSim = length(alpha_TrBB),seed = 2015)
save(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80,file="Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80.Rdata")


load("JASA Submission\\MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a.Rdata")

alpha_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a$target.parameters$alpha_sample),5000)
beta_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.a$target.parameters$beta_sample),5000)
cut_off <- 0.0025
hist(alpha_TrBB,breaks = 100)
hist(beta_TrBB,breaks = 100)

summary(alpha_TrBB)
summary(beta_TrBB)

set.seed(2015)
Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 10000,B = 8000,M = 40,m = 5,b = 13,sensitivity_values = c(0.7,0.8,0.9,0.95,1.0),NSim = length(alpha_TrBB),seed = 2015)
save(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80,file="Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80.Rdata")



load("JASA Submission\\MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b.Rdata")

alpha_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b$target.parameters$alpha_sample),5000)
beta_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.b$target.parameters$beta_sample),5000)
cut_off <- 0.005
hist(alpha_TrBB,breaks = 100)
hist(beta_TrBB,breaks = 100)

summary(alpha_TrBB)
summary(beta_TrBB)

set.seed(2015)
Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 10000,B = 8000,M = 40,m = 5,b = 13,sensitivity_values = c(0.7,0.8,0.9,0.95,1.0),NSim = length(alpha_TrBB),seed = 2015)
save(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80,file="Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80.Rdata")


load("JASA Submission\\MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c.Rdata")

alpha_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c$target.parameters$alpha_sample),5000)
beta_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c$target.parameters$beta_sample),5000)
cut_off <- 0.01
hist(alpha_TrBB,breaks = 100)
hist(beta_TrBB,breaks = 100)

summary(alpha_TrBB)
summary(beta_TrBB)

set.seed(2015)
Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 10000,B = 8000,M = 40,m = 5,b = 13,sensitivity_values = c(0.7,0.8,0.9,0.95,1.0),NSim = length(alpha_TrBB),seed = 2015)
save(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80,file="Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80.Rdata")



load("JASA Submission\\MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d.Rdata")

alpha_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d$target.parameters$alpha_sample),5000)
beta_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.d$target.parameters$beta_sample),5000)
cut_off <- 0.02
hist(alpha_TrBB,breaks = 100)
hist(beta_TrBB,breaks = 100)

summary(alpha_TrBB)
summary(beta_TrBB)

set.seed(2015)
Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 10000,B = 8000,M = 40,m = 5,b = 13,sensitivity_values = c(0.7,0.8,0.9,0.95,1.0),NSim = length(alpha_TrBB),seed = 2015)
save(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80,file="Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80.Rdata")


# Load the Model-based Simulation results ----------------------------
# The Simulated Data are not saved in GitHub due to restriction of data size
# Extract and compute proportions of simulated infection outcomes
# across different assumed k values (Sn = 0, 0.0025, 0.005, 0.01, 0.02).
# across different assumed sensitivity values (Sn = 0.7, 0.8, 0.9, 0.95, 1.0).
# The results provide the relative frequency distribution of
# infection types under each scenario.
# -------------------------------------------------------------------#

load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80.Rdata")
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80$`0.7`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80$`0.8`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80$`1`$Simulation_Infection_Ty_Results))

load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80.Rdata")
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80$`0.7`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80$`0.8`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80$`1`$Simulation_Infection_Ty_Results))

load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80.Rdata")
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80$`0.7`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80$`0.8`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80$`1`$Simulation_Infection_Ty_Results))

load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80.Rdata")
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`0.7`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`0.8`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`1`$Simulation_Infection_Ty_Results))

load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80.Rdata")
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80$`0.7`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80$`0.8`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80$`1`$Simulation_Infection_Ty_Results))


# Summary Statistics --------------------#
Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80 <- NULL
Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80 <- NULL
Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80 <- NULL
Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80 <- NULL
Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80 <- NULL

prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.7`$Simulation_Infection_Ty_Results)) # TP, TN and FN : 0.01204060 0.94853272 0.03942668
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.8`$Simulation_Infection_Ty_Results)) # TP, TN and FN : 0.01263662 0.94853272 0.03883066
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.1`$Simulation_Infection_Ty_Results))

# TP, TN and FN : 0.01269084 0.94981886 0.03749030
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect$`0.7`$Simulation_Infection_Ty_Results))
# TP, TN and FN : 0.01267208 0.97964040 0.00768752
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect$`0.7`$Simulation_Infection_Ty_Results))
# TP, TN and FN : 0.01268842 0.98199412 0.00531746
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.7`$Simulation_Infection_Ty_Results))
# TP, TN and FN : 0.01273912 0.98397560 0.00328528
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.7`$Simulation_Infection_Ty_Results))
# TP, TN and FN : 0.01297678 0.98530724 0.00171598


# -------------------------------------------------------------------
# Figures for Community Risk - Main Manuscript
# -------------------------------------------------------------------
# Simulation parameters are defined for generating figures on
# community-level infection risk under imperfect testing conditions.
#
# Key specifications:
#   - Delta = 0.80 (represents imperfect diagnostic sensitivity)
#   - k     = 0.00, 0.01, 0.02 (levels of overdispersion considered)
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Simulation Results Processing for Community Risk Figures
# -------------------------------------------------------------------
# This script iterates over selected values of:
#   - m (number of pools),
#   - b (pool size),
#   - delta (test sensitivity parameter),
#   - cutoff_values (overdispersion parameter settings).
#
# For each parameter combination, the corresponding simulation dataset
# is identified, validated, and processed. Results are stored in a
# structured list with descriptive keys for subsequent analysis.
# -------------------------------------------------------------------

# Function: process_tx_results --------------------------------------#
#
#
# PURPOSE:
#   Process simulation results for *Tx* (true infected units) and corresponding
#   inspection outcomes, given a specific sensitivity index.
#
# PROCESS FLOW:
#   1. Extract Tx (true infections) and inspection results from the simulation.
#   2. Convert Tx to numeric, and record the sensitivity used.
#   3. Convert inspection outcomes into labeled categories (TP, TN, FN).
#   4. Create Tx_range → classify Tx into meaningful ranges (0, 1–100, 101–500, …).
#   5. Convert Tx_range into an ordered factor for consistent plotting/analysis.
#   6. Return the processed dataframe.
#
# OUTPUT:
#   A dataframe with columns:
#     - Tx (numeric true infection counts)
#     - Inspection.Results (factor: TP/TN/FN)
#     - Sensitivity (which sensitivity setting was applied)
#     - Tx_range (categorical range of Tx values)
# -----------------------------------------------------------------------------

process_tx_results <- function(sim_data, index) {

  # ---- Step 1: Extract Tx and Inspection Results from the simulation ----
  df <- as.data.frame(cbind(
    Tx = as.vector(sim_data[[as.character(index)]]$Simulation_Tx),
    Inspection.Results = as.vector(sim_data[[as.character(index)]]$Simulation_Infection_Tx_Results)
  ))

  # ---- Step 2: Convert Tx to numeric and attach Sensitivity ----
  df$Tx <- as.numeric(paste(df$Tx))
  df$Sensitivity <- as.character(index)

  # ---- Step 3: Relabel Inspection Results ----
  # Convert to factor with descriptive labels: TP = True Positive, TN = True Negative, FN = False Negative
  df$Inspection.Results <- factor(
    as.factor(sim_data[[as.character(index)]]$Simulation_Infection_Ty_Results),
    labels = c("TP", "TN", "FN")
  )

  # ---- Step 4: Create Tx_range column (categorical bins for Tx values) ----
  df$Tx_range <- "0"
  df$Tx_range[df$Tx <= 100 & df$Tx > 0]    <- "1-100"
  df$Tx_range[df$Tx <= 500 & df$Tx > 100]  <- "101-500"
  df$Tx_range[df$Tx <= 1000 & df$Tx > 500] <- "501-1000"
  df$Tx_range[df$Tx <= 2000 & df$Tx > 1000] <- "1001-2000"
  df$Tx_range[df$Tx <= 3000 & df$Tx > 2000] <- "2001-3000"
  df$Tx_range[df$Tx <= 4000 & df$Tx > 3000] <- "3001-4000"
  df$Tx_range[df$Tx <= 5000 & df$Tx > 4000] <- "4001-5000"
  df$Tx_range[df$Tx <= 10000 & df$Tx > 5000] <- "5001-10000"
  df$Tx_range[df$Tx <= 15000 & df$Tx > 10001] <- "10001-15000"
  df$Tx_range[df$Tx <= 30000 & df$Tx > 15001] <- "15001-30000"
  df$Tx_range[df$Tx > 30000] <- "30000+"

  # ---- Step 5: Make Tx_range an ordered factor ----
  df$Tx_range <- factor(
    as.factor(df$Tx_range),
    levels = c("0", "1-100", "101-500", "501-1000", "1001-2000", "2001-3000",
               "3001-4000", "4001-5000", "5001-10000", "10001-15000", "15001-30000", "30000+")
  )

  # ---- Step 6: Return processed dataframe ----
  return(df)
}


# Define model parameters
m_values <- c(5)          # Number of pools
b_values <- c(13)         # Pool size
delta_values <- c(0.8)    # Imperfect sensitivity case only
cutoff_values <- c("_MH_00","_MH_0025","_MH_005","_MH_01","_MH_02")
# cutoff_values <- c("_MH_00","_MH_01","_MH_02")
# "_MH_00"   -> k = 0.00
# "_MH_0025" -> k = 0.0025
# "_MH_005"  -> k = 0.005
# "_MH_01"   -> k = 0.01
# "_MH_02"   -> k = 0.02

# Initialize an empty list to store processed results
results_list <- list()

# Iterate through all combinations of parameters
for (m in m_values) {
  for (b in b_values) {
    for (c in cutoff_values) {
      for (delta in delta_values) {

        # Dynamically construct the dataset name
        dataset_name <- paste0("Sim_Tx_TrBB_Prawn_m", m, "_b", b, c, "_Imperfect")

        # Verify dataset availability before processing
        if (exists(dataset_name)) {

          # Process simulation results using the custom function
          result <- process_tx_results(get(dataset_name), delta)

          # Store results in the list with a descriptive key
          key <- paste0("Sn_", delta, "_m_", m, "_b_", b, c)
          results_list[[key]] <- result

        } else {
          # Display a message if the dataset is missing
          message("Dataset ", dataset_name, " not found. Skipping...")
        }
      }
    }
  }
}



# -------------------------------------------------------------------
# Aggregating Processed Simulation Results into a Combined Data Frame
# -------------------------------------------------------------------
# This script:
#   1. Iterates through all combinations of m, b, cutoff, and sensitivity (Sn).
#   2. Extracts corresponding processed results from results_list.
#   3. Creates contingency tables for Tx range vs. Inspection Results.
#   4. Combines all results into a single standardized data frame.
#   5. Normalizes frequencies and computes percentages.
#   6. Saves the final combined dataset for subsequent analyses.
# -------------------------------------------------------------------

# Initialize an empty list to store intermediate data frames
df_list <- list()

# Nested loops over number of items (m), number of bags (b), cutoff values, and sensitivity (Sn)
for (m in m_values) {
  for (b in b_values) {
    for (c in cutoff_values) {
      for (Sn in delta_values) {

        # Construct the key used to locate the results in results_list
        df_key <- paste0("Sn_", Sn, "_m_", m, "_b_", b, c)

        # Check if this key exists in results_list
        if (df_key %in% names(results_list)) {

          # Define a descriptive name for the extracted dataframe
          df_result_name <- paste("df.perfect.Sn", Sn, "m", m, "b", b, c, sep = ".")

          # Extract the corresponding dataframe from results_list
          df_data <- results_list[[df_key]]

          # Create a contingency table of Tx_range vs. Inspection.Results
          df <- as.data.frame(table(df_data$Tx_range, df_data$Inspection.Results))
          names(df) <- c("Tx", "Outcome", "Freq")  # Rename columns

          # Add metadata columns for traceability
          df$m <- as.character(m)
          df$b <- as.character(b)
          df$Sn <- as.character(Sn)
          df$cut <- as.character(c)

          # Store the processed dataframe into df_list
          df_list[[df_result_name]] <- df

        } else {
          # Display a message if the key is missing
          message("Key ", df_key, " not found in results_list. Skipping...")
        }
      }
    }
  }
}

# -------------------------------------------------------------------
# Combine all data frames into one standardized data frame
# -------------------------------------------------------------------
df_combined <- do.call(rbind, lapply(names(df_list), function(name) {
  df <- df_list[[name]]
  df$source <- name  # Add a column to record the source dataframe name
  return(df)
}))

# -------------------------------------------------------------------
# Post-processing: add derived variables and normalize frequencies
# -------------------------------------------------------------------

# Label counts as "zero" or "nonzero" based on Tx value
df_combined$Counts <- ifelse(df_combined$Tx == "0", "zero", "nonzero")

# Count unique values for scaling frequencies
m_counts <- length(unique(df_combined$m))
b_counts <- length(unique(df_combined$b))
sn_counts <- length(unique(df_combined$Sn))
cutoff_counts <- length(unique(df_combined$cut))

# Normalize frequencies to ensure comparability across scenarios
df_combined$Freq <- (df_combined$Freq / sum(df_combined$Freq)) *
  (1000 * m_counts * b_counts * sn_counts * cutoff_counts)

# Compute percentages based on normalized frequencies
df_combined$Perc <- (df_combined$Freq / 1000) * 100

# -------------------------------------------------------------------
# Save the final combined dataset for future analyses
# -------------------------------------------------------------------
save(df_combined, file = "df_combined_Model_Based_Simulation_Sn.0.80_varying_k.Rdata")


# Figures for Community Risk - Main Manuscript and Supplementary Material --------------

library(ggplot2)
library(ggpubr)

load("JASA Submission\\df_combined_Model_Based_Simulation_Sn.0.80_varying_k.Rdata")

m_values <- c(5)
b_values <- c(13)

# Initialize an empty list to store the plots
plot_list_imperfect <- list()
cutoff_values_trBB <- c("_MH_00","_MH_0025","_MH_005","_MH_01","_MH_02") # DAFF Presentation
# cutoff_values_trBB <- c("_MH_00","_MH_01","_MH_02") # DAFF Presentation
# Loop through all combinations of m and b
for (m in m_values) {
  for (b in b_values) {
    for (c in cutoff_values_trBB) {

      # Map the cutoff value to its corresponding label
      cutoff_label <- switch(c,
                             "_MH_00" = "0.00",
                             "_MH_0025" = "0.0025",
                             "_MH_005" = "0.005",
                             "_MH_01" = "0.01",
                             "_MH_02" = "0.02",
                             "_Unknown")   # Default label for unexpected values
      # Filter the data based on current m and b values
      df_nonzero <- df_combined[df_combined$Counts == "nonzero" & df_combined$m == m & df_combined$b == b & df_combined$cut == c, ]
      df_zero <- df_combined[df_combined$Counts == "zero" & df_combined$m == m & df_combined$b == b & df_combined$cut == c, ]

      # max_freq <- ceiling(max(df_nonzero$Freq))+0.05
      max_freq <- 23+0.05

      # Plot for nonzero counts
      nonzero <- ggplot(df_nonzero, aes(x = Sn, y = Freq, fill = Outcome)) +
        geom_col(position = "stack") +
        facet_wrap(~ Tx, nrow = 1) +
        scale_fill_manual(values = c("TP" = "blue", "TN" = "green", "FN" = "red")) +
        labs(x = NULL, title = paste("Tx > 0")) +
        scale_y_continuous(
          name = "Frequency (in 1000s)",
          breaks = seq(0, max_freq, by = 2),limits = c(0,max_freq)
        )+
        theme(
          axis.text.x = element_blank(),   # Removes x-axis tick labels
          axis.ticks.x = element_blank(),  # Removes x-axis tick marks (optional)
          strip.text = element_text(size = 8,face = "bold")  # Reduce facet label font size to 8
        )

      # Plot for zero counts
      zero <- ggplot(df_zero, aes(x = Sn, y = Freq, fill = Outcome)) +
        geom_col(position = "stack") +
        scale_fill_manual(values = c("TP" = "blue", "TN" = "green", "FN" = "red")) +
        facet_wrap(~ Tx, nrow = 1) +
        labs(x = NULL, title = paste("Tx = 0")) +
        scale_y_continuous(
          name = "Frequency (in 1000s)",
          breaks = seq(0, 1000, by = 100),limits = c(0,1000)
        )+
        theme(
          axis.text.x = element_blank(),   # Removes x-axis tick labels
          axis.ticks.x = element_blank(),  # Removes x-axis tick marks (optional)
          strip.text = element_text(size = 8,face = "bold")  # Reduce facet label font size to 8
        )

      # Combine the zero and nonzero plots

      plot_m_b <- annotate_figure(
        ggarrange(zero, nonzero, nrow = 1, common.legend = TRUE, legend = "top",
                  widths = c(0.15, 0.85)),
        top = text_grob(
          # bquote("Number of batches by outcome & infected prawns (" * T[X] * "): k = " * .(cutoff_label)),  # Use bquote() for math formatting
          bquote("Minimum Positive Propensity (" * italic(k) * ") : " * .(cutoff_label)),  # Use bquote() for math formatting
          face = "bold",
          size = 14
        )
      )

      # Store the plot in the list with a unique name
      plot_name <- paste("plot_m", m, "b", b, c, sep = "_")
      plot_list_imperfect[[plot_name]] <- plot_m_b
    }
  }
}



ggarrange(plot_list_imperfect$plot_m_5_b_13__MH_00,
          plot_list_imperfect$plot_m_5_b_13__MH_01,
          plot_list_imperfect$plot_m_5_b_13__MH_02,
          nrow=3,ncol=1,common.legend = TRUE)

ggsave("Test_Outcome_m5_b13_Imperfect_Manuscript.png",width = 10,height = 10,units = "in")

ggarrange(plot_list_imperfect$plot_m_5_b_13__MH_00,
          plot_list_imperfect$plot_m_5_b_13__MH_0025,
          plot_list_imperfect$plot_m_5_b_13__MH_005,
          plot_list_imperfect$plot_m_5_b_13__MH_01,
          plot_list_imperfect$plot_m_5_b_13__MH_02,
          nrow=5,ncol=1,common.legend = TRUE)

ggsave("Test_Outcome_m5_b13_Imperfect_Manuscript_SupMat.png",width = 10,height = 16,units = "in")



# Supplementary Table for Community Risk - Supplementary Material --------------

load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80 <- NULL

load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect_Sn_80 <- NULL
load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect_Sn_80 <- NULL
load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80 <- NULL
load("JASA Submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80 <- NULL


mean_LKG_B_00_perf <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`1`)
mean_LKG_B_00_imperf_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.95`)
mean_LKG_B_00_imperf_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.9`)
mean_LKG_B_00_imperf_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.8`)
mean_LKG_B_00_imperf_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.7`)

mean_LKG_B_0025_perf <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect$`1`)
mean_LKG_B_0025_imperf_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect$`0.95`)
mean_LKG_B_0025_imperf_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect$`0.9`)
mean_LKG_B_0025_imperf_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect$`0.8`)
mean_LKG_B_0025_imperf_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_0025_Imperfect$`0.7`)

mean_LKG_B_005_perf <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect$`1`)
mean_LKG_B_005_imperf_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect$`0.95`)
mean_LKG_B_005_imperf_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect$`0.9`)
mean_LKG_B_005_imperf_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect$`0.8`)
mean_LKG_B_005_imperf_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_005_Imperfect$`0.7`)

mean_LKG_B_01_perf <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`1`)
mean_LKG_B_01_imperf_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.95`)
mean_LKG_B_01_imperf_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.9`)
mean_LKG_B_01_imperf_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.8`)
mean_LKG_B_01_imperf_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.7`)

mean_LKG_B_02_perf <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`1`)
mean_LKG_B_02_imperf_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.95`)
mean_LKG_B_02_imperf_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.9`)
mean_LKG_B_02_imperf_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.8`)
mean_LKG_B_02_imperf_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.7`)


mean_LKG_df_perf_1 <- data.frame(
  E.Li = c(
    mean_LKG_B_00_perf$E.Li,
    mean_LKG_B_0025_perf$E.Li,
    mean_LKG_B_005_perf$E.Li,
    mean_LKG_B_01_perf$E.Li,
    mean_LKG_B_02_perf$E.Li
  ),
  P.Li = c(
    mean_LKG_B_00_perf$P.Li,
    mean_LKG_B_0025_perf$P.Li,
    mean_LKG_B_005_perf$P.Li,
    mean_LKG_B_01_perf$P.Li,
    mean_LKG_B_02_perf$P.Li
  ),
  Median.Li = c(
    mean_LKG_B_00_perf$Median.Li,
    mean_LKG_B_0025_perf$Median.Li,
    mean_LKG_B_005_perf$Median.Li,
    mean_LKG_B_01_perf$Median.Li,
    mean_LKG_B_02_perf$Median.Li
  ),
  Pseudo.E.Li.median = c(
    mean_LKG_B_00_perf$Pseudo.E.Li.median,
    mean_LKG_B_0025_perf$Pseudo.E.Li.median,
    mean_LKG_B_005_perf$Pseudo.E.Li.median,
    mean_LKG_B_01_perf$Pseudo.E.Li.median,
    mean_LKG_B_02_perf$Pseudo.E.Li.median
  ),
  cutoff = factor(rep(c(0.00, 0.0025, 0.005, 0.01, 0.02), each = length(mean_LKG_B_0025_perf$P.Li))),
  sensitivity = rep(c("1.0"), times = 5*length(mean_LKG_B_0025_perf$P.Li))
)


P.Li.perf <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_perf_1, FUN = mean)
P.Li.perf.CI <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_perf_1, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E.Li.perf <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_perf_1, FUN = mean)
E.Li.perf.CI <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_perf_1, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Median.Li.perf <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_perf_1, FUN = median)
Median.Li.perf.CI <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_perf_1, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Pseudo.E.Li.median.perf <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_perf_1, FUN = median)
Pseudo.E.Li.median.perf.CI <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_perf_1, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P.E.Li.perf.sen <- as.data.frame(E.Li.perf)
P.E.Li.perf.sen$E.Li.ll <- E.Li.perf.CI$E.Li[,1]
P.E.Li.perf.sen$E.Li.ul <- E.Li.perf.CI$E.Li[,2]
P.E.Li.perf.sen$P.Li <- P.Li.perf$P.Li
P.E.Li.perf.sen$P.Li.ll <- P.Li.perf.CI$P.Li[,1]
P.E.Li.perf.sen$P.Li.ul <- P.Li.perf.CI$P.Li[,2]
P.E.Li.perf.sen$Median.Li <- Median.Li.perf$Median.Li
P.E.Li.perf.sen$Median.Li.ll <- Median.Li.perf.CI$Median.Li[,1]
P.E.Li.perf.sen$Median.Li.ul <- Median.Li.perf.CI$Median.Li[,2]
P.E.Li.perf.sen$Pseudo.E.Li.median <- Pseudo.E.Li.median.perf$Pseudo.E.Li.median
P.E.Li.perf.sen$Pseudo.E.Li.median.ll <- Pseudo.E.Li.median.perf.CI$Pseudo.E.Li.median[,1]
P.E.Li.perf.sen$Pseudo.E.Li.median.ul <- Pseudo.E.Li.median.perf.CI$Pseudo.E.Li.median[,2]


mean_LKG_df_imperf_0.7 <- data.frame(
  E.Li = c(
    mean_LKG_B_00_imperf_0.7$E.Li,
    mean_LKG_B_0025_imperf_0.7$E.Li,
    mean_LKG_B_005_imperf_0.7$E.Li,
    mean_LKG_B_01_imperf_0.7$E.Li,
    mean_LKG_B_02_imperf_0.7$E.Li
  ),
  P.Li = c(
    mean_LKG_B_00_imperf_0.7$P.Li,
    mean_LKG_B_0025_imperf_0.7$P.Li,
    mean_LKG_B_005_imperf_0.7$P.Li,
    mean_LKG_B_01_imperf_0.7$P.Li,
    mean_LKG_B_02_imperf_0.7$P.Li
  ),
  Median.Li = c(
    mean_LKG_B_00_imperf_0.7$Median.Li,
    mean_LKG_B_0025_imperf_0.7$Median.Li,
    mean_LKG_B_005_imperf_0.7$Median.Li,
    mean_LKG_B_01_imperf_0.7$Median.Li,
    mean_LKG_B_02_imperf_0.7$Median.Li
  ),
  Pseudo.E.Li.median = c(
    mean_LKG_B_00_imperf_0.7$Pseudo.E.Li.median,
    mean_LKG_B_0025_imperf_0.7$Pseudo.E.Li.median,
    mean_LKG_B_005_imperf_0.7$Pseudo.E.Li.median,
    mean_LKG_B_01_imperf_0.7$Pseudo.E.Li.median,
    mean_LKG_B_02_imperf_0.7$Pseudo.E.Li.median
  ),
  cutoff = factor(rep(c(0.00, 0.0025, 0.005, 0.01, 0.02), each = length(mean_LKG_B_0025_imperf_0.7$P.Li))),
  sensitivity = rep(c("0.7"), times = 5*length(mean_LKG_B_0025_imperf_0.7$P.Li))
)

P.Li.imperf.70 <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.7, FUN = mean)
P.Li.imperf.CI.70 <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.7, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E.Li.imperf.70 <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.7, FUN = mean)
E.Li.imperf.CI.70 <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.7, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Median.Li.imperf.70 <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.7, FUN = median)
Median.Li.imperf.CI.70 <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.7, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Pseudo.E.Li.median.imperf.70 <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.7, FUN = median)
Pseudo.E.Li.median.imperf.CI.70 <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.7, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P.E.Li.imperf.sen.70 <- as.data.frame(E.Li.imperf.70)
P.E.Li.imperf.sen.70$E.Li.ll <- E.Li.imperf.CI.70$E.Li[,1]
P.E.Li.imperf.sen.70$E.Li.ul <- E.Li.imperf.CI.70$E.Li[,2]
P.E.Li.imperf.sen.70$P.Li <- P.Li.imperf.70$P.Li
P.E.Li.imperf.sen.70$P.Li.ll <- P.Li.imperf.CI.70$P.Li[,1]
P.E.Li.imperf.sen.70$P.Li.ul <- P.Li.imperf.CI.70$P.Li[,2]
P.E.Li.imperf.sen.70$Median.Li <- Median.Li.imperf.70$Median.Li
P.E.Li.imperf.sen.70$Median.Li.ll <- Median.Li.imperf.CI.70$Median.Li[,1]
P.E.Li.imperf.sen.70$Median.Li.ul <- Median.Li.imperf.CI.70$Median.Li[,2]
P.E.Li.imperf.sen.70$Pseudo.E.Li.median <- Pseudo.E.Li.median.imperf.70$Pseudo.E.Li.median
P.E.Li.imperf.sen.70$Pseudo.E.Li.median.ll <- Pseudo.E.Li.median.imperf.CI.70$Pseudo.E.Li.median[,1]
P.E.Li.imperf.sen.70$Pseudo.E.Li.median.ul <- Pseudo.E.Li.median.imperf.CI.70$Pseudo.E.Li.median[,2]



mean_LKG_df_imperf_0.8 <- data.frame(
  E.Li = c(
    mean_LKG_B_00_imperf_0.8$E.Li,
    mean_LKG_B_0025_imperf_0.8$E.Li,
    mean_LKG_B_005_imperf_0.8$E.Li,
    mean_LKG_B_01_imperf_0.8$E.Li,
    mean_LKG_B_02_imperf_0.8$E.Li
  ),
  P.Li = c(
    mean_LKG_B_00_imperf_0.8$P.Li,
    mean_LKG_B_0025_imperf_0.8$P.Li,
    mean_LKG_B_005_imperf_0.8$P.Li,
    mean_LKG_B_01_imperf_0.8$P.Li,
    mean_LKG_B_02_imperf_0.8$P.Li
  ),
  Median.Li = c(
    mean_LKG_B_00_imperf_0.8$Median.Li,
    mean_LKG_B_0025_imperf_0.8$Median.Li,
    mean_LKG_B_005_imperf_0.8$Median.Li,
    mean_LKG_B_01_imperf_0.8$Median.Li,
    mean_LKG_B_02_imperf_0.8$Median.Li
  ),
  Pseudo.E.Li.median = c(
    mean_LKG_B_00_imperf_0.8$Pseudo.E.Li.median,
    mean_LKG_B_0025_imperf_0.8$Pseudo.E.Li.median,
    mean_LKG_B_005_imperf_0.8$Pseudo.E.Li.median,
    mean_LKG_B_01_imperf_0.8$Pseudo.E.Li.median,
    mean_LKG_B_02_imperf_0.8$Pseudo.E.Li.median
  ),
  cutoff = factor(rep(c(0.00, 0.0025, 0.005, 0.01, 0.02), each = length(mean_LKG_B_0025_imperf_0.8$P.Li))),
  sensitivity = rep(c("0.8"), times = 5*length(mean_LKG_B_0025_imperf_0.8$P.Li))
)

P.Li.imperf.80 <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.8, FUN = mean)
P.Li.imperf.CI.80 <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.8, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E.Li.imperf.80 <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.8, FUN = mean)
E.Li.imperf.CI.80 <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.8, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Median.Li.imperf.80 <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.8, FUN = median)
Median.Li.imperf.CI.80 <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.8, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Pseudo.E.Li.median.imperf.80 <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.8, FUN = median)
Pseudo.E.Li.median.imperf.CI.80 <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.8, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P.E.Li.imperf.sen.80 <- as.data.frame(E.Li.imperf.80)
P.E.Li.imperf.sen.80$E.Li.ll <- E.Li.imperf.CI.80$E.Li[,1]
P.E.Li.imperf.sen.80$E.Li.ul <- E.Li.imperf.CI.80$E.Li[,2]
P.E.Li.imperf.sen.80$P.Li <- P.Li.imperf.80$P.Li
P.E.Li.imperf.sen.80$P.Li.ll <- P.Li.imperf.CI.80$P.Li[,1]
P.E.Li.imperf.sen.80$P.Li.ul <- P.Li.imperf.CI.80$P.Li[,2]
P.E.Li.imperf.sen.80$Median.Li <- Median.Li.imperf.80$Median.Li
P.E.Li.imperf.sen.80$Median.Li.ll <- Median.Li.imperf.CI.80$Median.Li[,1]
P.E.Li.imperf.sen.80$Median.Li.ul <- Median.Li.imperf.CI.80$Median.Li[,2]
P.E.Li.imperf.sen.80$Pseudo.E.Li.median <- Pseudo.E.Li.median.imperf.80$Pseudo.E.Li.median
P.E.Li.imperf.sen.80$Pseudo.E.Li.median.ll <- Pseudo.E.Li.median.imperf.CI.80$Pseudo.E.Li.median[,1]
P.E.Li.imperf.sen.80$Pseudo.E.Li.median.ul <- Pseudo.E.Li.median.imperf.CI.80$Pseudo.E.Li.median[,2]



mean_LKG_df_imperf_0.9 <- data.frame(
  E.Li = c(
    mean_LKG_B_00_imperf_0.9$E.Li,
    mean_LKG_B_0025_imperf_0.9$E.Li,
    mean_LKG_B_005_imperf_0.9$E.Li,
    mean_LKG_B_01_imperf_0.9$E.Li,
    mean_LKG_B_02_imperf_0.9$E.Li
  ),
  P.Li = c(
    mean_LKG_B_00_imperf_0.9$P.Li,
    mean_LKG_B_0025_imperf_0.9$P.Li,
    mean_LKG_B_005_imperf_0.9$P.Li,
    mean_LKG_B_01_imperf_0.9$P.Li,
    mean_LKG_B_02_imperf_0.9$P.Li
  ),
  Median.Li = c(
    mean_LKG_B_00_imperf_0.9$Median.Li,
    mean_LKG_B_0025_imperf_0.9$Median.Li,
    mean_LKG_B_005_imperf_0.9$Median.Li,
    mean_LKG_B_01_imperf_0.9$Median.Li,
    mean_LKG_B_02_imperf_0.9$Median.Li
  ),
  Pseudo.E.Li.median = c(
    mean_LKG_B_00_imperf_0.9$Pseudo.E.Li.median,
    mean_LKG_B_0025_imperf_0.9$Pseudo.E.Li.median,
    mean_LKG_B_005_imperf_0.9$Pseudo.E.Li.median,
    mean_LKG_B_01_imperf_0.9$Pseudo.E.Li.median,
    mean_LKG_B_02_imperf_0.9$Pseudo.E.Li.median
  ),
  cutoff = factor(rep(c(0.00, 0.0025, 0.005, 0.01, 0.02), each = length(mean_LKG_B_0025_imperf_0.9$P.Li))),
  sensitivity = rep(c("0.9"), times = 5*length(mean_LKG_B_0025_imperf_0.9$P.Li))
)

P.Li.imperf.90 <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.9, FUN = mean)
P.Li.imperf.CI.90 <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.9, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E.Li.imperf.90 <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.9, FUN = mean)
E.Li.imperf.CI.90 <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.9, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Median.Li.imperf.90 <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.9, FUN = median)
Median.Li.imperf.CI.90 <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.9, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Pseudo.E.Li.median.imperf.90 <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.9, FUN = median)
Pseudo.E.Li.median.imperf.CI.90 <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.9, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P.E.Li.imperf.sen.90 <- as.data.frame(E.Li.imperf.90)
P.E.Li.imperf.sen.90$E.Li.ll <- E.Li.imperf.CI.90$E.Li[,1]
P.E.Li.imperf.sen.90$E.Li.ul <- E.Li.imperf.CI.90$E.Li[,2]
P.E.Li.imperf.sen.90$P.Li <- P.Li.imperf.90$P.Li
P.E.Li.imperf.sen.90$P.Li.ll <- P.Li.imperf.CI.90$P.Li[,1]
P.E.Li.imperf.sen.90$P.Li.ul <- P.Li.imperf.CI.90$P.Li[,2]
P.E.Li.imperf.sen.90$Median.Li <- Median.Li.imperf.90$Median.Li
P.E.Li.imperf.sen.90$Median.Li.ll <- Median.Li.imperf.CI.90$Median.Li[,1]
P.E.Li.imperf.sen.90$Median.Li.ul <- Median.Li.imperf.CI.90$Median.Li[,2]
P.E.Li.imperf.sen.90$Pseudo.E.Li.median <- Pseudo.E.Li.median.imperf.90$Pseudo.E.Li.median
P.E.Li.imperf.sen.90$Pseudo.E.Li.median.ll <- Pseudo.E.Li.median.imperf.CI.90$Pseudo.E.Li.median[,1]
P.E.Li.imperf.sen.90$Pseudo.E.Li.median.ul <- Pseudo.E.Li.median.imperf.CI.90$Pseudo.E.Li.median[,2]



mean_LKG_df_imperf_0.95 <- data.frame(
  E.Li = c(
    mean_LKG_B_00_imperf_0.95$E.Li,
    mean_LKG_B_0025_imperf_0.95$E.Li,
    mean_LKG_B_005_imperf_0.95$E.Li,
    mean_LKG_B_01_imperf_0.95$E.Li,
    mean_LKG_B_02_imperf_0.95$E.Li
  ),
  P.Li = c(
    mean_LKG_B_00_imperf_0.95$P.Li,
    mean_LKG_B_0025_imperf_0.95$P.Li,
    mean_LKG_B_005_imperf_0.95$P.Li,
    mean_LKG_B_01_imperf_0.95$P.Li,
    mean_LKG_B_02_imperf_0.95$P.Li
  ),
  Median.Li = c(
    mean_LKG_B_00_imperf_0.95$Median.Li,
    mean_LKG_B_0025_imperf_0.95$Median.Li,
    mean_LKG_B_005_imperf_0.95$Median.Li,
    mean_LKG_B_01_imperf_0.95$Median.Li,
    mean_LKG_B_02_imperf_0.95$Median.Li
  ),
  Pseudo.E.Li.median = c(
    mean_LKG_B_00_imperf_0.95$Pseudo.E.Li.median,
    mean_LKG_B_0025_imperf_0.95$Pseudo.E.Li.median,
    mean_LKG_B_005_imperf_0.95$Pseudo.E.Li.median,
    mean_LKG_B_01_imperf_0.95$Pseudo.E.Li.median,
    mean_LKG_B_02_imperf_0.95$Pseudo.E.Li.median
  ),
  cutoff = factor(rep(c(0.00, 0.0025, 0.005, 0.01, 0.02), each = length(mean_LKG_B_0025_imperf_0.95$P.Li))),
  sensitivity = rep(c("0.95"), times = 5*length(mean_LKG_B_0025_imperf_0.95$P.Li))
)

P.Li.imperf.95 <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.95, FUN = mean)
P.Li.imperf.CI.95 <- aggregate(P.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.95, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E.Li.imperf.95 <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.95, FUN = mean)
E.Li.imperf.CI.95 <- aggregate(E.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.95, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Median.Li.imperf.95 <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.95, FUN = median)
Median.Li.imperf.CI.95 <- aggregate(Median.Li ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.95, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
Pseudo.E.Li.median.imperf.95 <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.95, FUN = median)
Pseudo.E.Li.median.imperf.CI.95 <- aggregate(Pseudo.E.Li.median ~ cutoff + sensitivity, data = mean_LKG_df_imperf_0.95, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P.E.Li.imperf.sen.95 <- as.data.frame(E.Li.imperf.95)
P.E.Li.imperf.sen.95$E.Li.ll <- E.Li.imperf.CI.95$E.Li[,1]
P.E.Li.imperf.sen.95$E.Li.ul <- E.Li.imperf.CI.95$E.Li[,2]
P.E.Li.imperf.sen.95$P.Li <- P.Li.imperf.95$P.Li
P.E.Li.imperf.sen.95$P.Li.ll <- P.Li.imperf.CI.95$P.Li[,1]
P.E.Li.imperf.sen.95$P.Li.ul <- P.Li.imperf.CI.95$P.Li[,2]
P.E.Li.imperf.sen.95$Median.Li <- Median.Li.imperf.95$Median.Li
P.E.Li.imperf.sen.95$Median.Li.ll <- Median.Li.imperf.CI.95$Median.Li[,1]
P.E.Li.imperf.sen.95$Median.Li.ul <- Median.Li.imperf.CI.95$Median.Li[,2]
P.E.Li.imperf.sen.95$Pseudo.E.Li.median <- Pseudo.E.Li.median.imperf.95$Pseudo.E.Li.median
P.E.Li.imperf.sen.95$Pseudo.E.Li.median.ll <- Pseudo.E.Li.median.imperf.CI.95$Pseudo.E.Li.median[,1]
P.E.Li.imperf.sen.95$Pseudo.E.Li.median.ul <- Pseudo.E.Li.median.imperf.CI.95$Pseudo.E.Li.median[,2]

E.P.Li.Community.Risk <- rbind(P.E.Li.perf.sen,P.E.Li.imperf.sen.95,P.E.Li.imperf.sen.90,P.E.Li.imperf.sen.80,P.E.Li.imperf.sen.70)

# write.csv(E.P.Li.Community.Risk,file="E.P.Li.Community.Risk.csv")
write.csv(E.P.Li.Community.Risk,file="E.P.Li.Community.Risk.Sn.80.csv")



# --------------------------------------------
# Create Supp Mat Table S4
# --------------------------------------------

E.P.Li.Community.Risk <- read.csv("E.P.Li.Community.Risk.Sn.80.csv",header = TRUE)

E.P.Li.Community.Risk <- E.P.Li.Community.Risk %>%
  mutate(across(
    .cols = where(is.numeric),
    .fns = ~ round(.x, 4)
  ))
E.P.Li.Community.Risk <- E.P.Li.Community.Risk %>%
  mutate(across(
    .cols = where(is.numeric) & !matches("^P\\.Li(\\.ll|\\.ul)?$"),
    .fns = ~ round(.x, 2)
  ))

digits_vec <- c(
  0,    # row names
  0,    # cutoff (factor — no decimals)
  0,    # inspection (character — no decimals)
  2,    # E.Li
  2,    # E.Li.ll
  2,    # E.Li.ul
  4,    # P.Li
  4,    # P.Li.ll
  4,    # P.Li.ul
  2,    # Median.Li
  2,    # Median.Li.ll
  2,    # Median.Li.ul
  2,    # Pseudo.E.Li.median
  2,    # Pseudo.E.Li.median.ll
  2     # Pseudo.E.Li.median.ul
)

print(xtable(E.P.Li.Community.Risk,digits = digits_vec),include.rownames=FALSE)


# Final Table: Only for Sensitivty:0.80

E.P.Li.Community.Risk.sn.80 <- E.P.Li.Community.Risk[E.P.Li.Community.Risk$sensitivity=="0.8",]

print(xtable(E.P.Li.Community.Risk.sn.80,digits = digits_vec),include.rownames=FALSE)

