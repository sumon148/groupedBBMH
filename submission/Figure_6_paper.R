# --------------------------------------------------------------------------- #
# Load Required Libraries
# --------------------------------------------------------------------------- #
library(groupedBBMH)  # Custom package for grouped beta-binomial MH modeling
library(ggplot2)
library(ggpubr)
library(xtable)

#-------------------------------------------------------------@
# Figure 6: JASA Paper
# Supplementary Table S5
# Community Risk: Vary number of sampled bags
# Model-based Simulation using MH-based TrBB model under Imperfect Testing
# Using Frozen-food imported data
# Exact BB Model: Bayesian method MH Algorithm
# Target Parameters: TP, TN, FN
# FN: Leakage
# P(Li>0): FN Batches
# E(Li): Number of contaminated items in FN Batches
# Delta=0.8, k=0.01, m=5, b=13, 20, 26
#-------------------------------------------------------------@

# -------------------------------------------------------------------
# NOTE: Do Not Run
# Model-based Simulation: MH-based TrBB model under Imperfect Test
# -------------------------------------------------------------------
# Purpose:
#   These lines reference pre-computed simulation results and should
#   not be executed directly during routine runs. They serve as a
#   record of the simulation setup and parameter configuration.
#
# Configuration:
#   - Delta (test sensitivity adjustment) = 0.80
#   - Lambda = 1
#   - k (contamination parameter) = 0.01
#   - b (number of bags) = 13, 20, 26
#
# Parameter Basis:
#   Beta distribution parameters used in this model were derived
#   under the specific assumption of delta = 0.80 and k = 0.01.
#
# ------------------------------------------------------------------- #

load("submission\\MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c.Rdata")

alpha_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c$target.parameters$alpha_sample),5000)
beta_TrBB <- sample(unlist(MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.c$target.parameters$beta_sample),5000)
cut_off <- 0.01
hist(alpha_TrBB,breaks = 100)
hist(beta_TrBB,breaks = 100)

summary(alpha_TrBB)
summary(beta_TrBB)

# b=13
set.seed(2015)
Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 10000,B = 8000,M = 40,m = 5,b = 13,sensitivity_values = c(0.7,0.8,0.9,0.95,1.0),NSim = length(alpha_TrBB),seed = 2015)
save(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80,file="Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80.Rdata")

# b=20
set.seed(2015)
Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 10000,B = 8000,M = 40,m = 5,b = 20,sensitivity_values = c(0.7,0.8,0.9,0.95,1.0),NSim = length(alpha_TrBB),seed = 2015)
save(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80,file="Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80.Rdata")

# b=26
set.seed(2015)
Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80 <- run_infection_simulations(alpha=alpha_TrBB,beta=beta_TrBB, cutoff = cut_off, D = 10000,B = 8000,M = 40,m = 5,b = 26,sensitivity_values = c(0.7,0.8,0.9,0.95,1.0),NSim = length(alpha_TrBB),seed = 2015)
save(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80,file="Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80.Rdata")


prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`0.7`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`0.8`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80$`1`$Simulation_Infection_Ty_Results))

prop.table(table(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80$`0.7`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80$`0.8`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80$`1`$Simulation_Infection_Ty_Results))

prop.table(table(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80$`0.7`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80$`0.8`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80$`0.9`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80$`0.95`$Simulation_Infection_Ty_Results))
prop.table(table(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80$`1`$Simulation_Infection_Ty_Results))


# Load the Model-based simulation -------------------------------------------
# Define parameter values for Figure 2: Community Risk under varying b=13,20,26 values
# Using model m=5 with Metropolis-Hastings (MH), imperfect sensitivity (Sn = 80%)

load("submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80.Rdata")
load("submission\\Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80.Rdata")
load("submission\\Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80.Rdata")

Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80


# Create Dataframe for creating plot --------------------------------------------------
# Define parameter values --------------------------------------------------#
# m_values: number of items inspected per bag (fixed at 5)
m_values <- c(5)

# b_values: different bag sizes used in the simulation
b_values <- c(13,20,26)

# delta_values: test sensitivity values
# Only imperfect sampling (Sn = 0.8) is considered in this run
# delta_values <- c(1.0, 0.7) # Uncomment for including perfect/imperfect
delta_values <- c(0.8)

# cutoff_values: suffix used in the dataset naming to indicate a specific model configuration
cutoff_values <- c("_MH_01")

# Initialize an empty list to store processed simulation results
results_list <- list()

# Loop over combinations of m, b, cutoff (c), and delta (Sn)
for (m in m_values) {
  for (b in b_values) {
    for (c in cutoff_values) {
      for (delta in delta_values) {
        # Construct dataset name string dynamically based on parameter values
        dataset_name <- paste0("Sim_Tx_TrBB_Prawn_m", m, "_b", b, c, "_Imperfect")

        # Check if the dataset exists in the current environment
        if (exists(dataset_name)) {
          # Apply processing function to simulation dataset
          result <- process_tx_results(get(dataset_name), delta)

          # Create a descriptive key for storing results
          key <- paste0("Sn_", delta, "_m_", m, "_b_", b, c)
          results_list[[key]] <- result
        } else {
          # Print a message if the expected dataset was not found
          message("Dataset ", dataset_name, " not found. Skipping...")
        }
      }
    }
  }
}

# -------------------------------------------------------------------------#
# Prepare output data frames for contingency table visualization ----------#

# Initialize empty list to store output data frames
df_list <- list()

# Loop again through the same parameter combinations
for (m in m_values) {
  for (b in b_values) {
    for (c in cutoff_values) {
      for (Sn in delta_values) {
        # Construct key to retrieve result from results_list
        df_key <- paste0("Sn_", Sn, "_m_", m, "_b_", b, c)

        # Check if the result exists
        if (df_key %in% names(results_list)) {
          # Define a name for the resulting data frame
          df_result_name <- paste("df.perfect.Sn", Sn, "m", m, "b", b, c, sep = ".")

          # Retrieve result data
          df_data <- results_list[[df_key]]

          # Generate contingency table from Tx status and inspection result
          df <- as.data.frame(table(df_data$Tx_range, df_data$Inspection.Results))
          names(df) <- c("Tx", "Outcome", "Freq")

          # Add parameter metadata to the table
          df$m <- as.character(m)
          df$b <- as.character(b)
          df$Sn <- as.character(Sn)
          df$cut <- as.character(c)

          # Store the processed table in list
          df_list[[df_result_name]] <- df
        } else {
          # If key is not found, skip and notify
          message("Key ", df_key, " not found in results_list. Skipping...")
        }
      }
    }
  }
}

# -------------------------------------------------------------------------#
# Combine all contingency tables into a single data frame -----------------#

df_combined <- do.call(rbind, lapply(names(df_list), function(name) {
  df <- df_list[[name]]
  df$source <- name  # Add column indicating source table name
  return(df)
}))

# Classify transmission status as zero or nonzero for easier grouping
df_combined$Counts <- ifelse(df_combined$Tx == "0", "zero", "nonzero")

# Count unique combinations of parameters for scaling purposes
m_counts <- length(unique(df_combined$m))
b_counts <- length(unique(df_combined$b))
sn_counts <- length(unique(df_combined$Sn))
cutoff_counts <- length(unique(df_combined$cut))

# Normalize frequencies across all parameter combinations
df_combined$Freq <- (df_combined$Freq / sum(df_combined$Freq)) *
  (1000 * m_counts * b_counts * sn_counts * cutoff_counts)

# Convert frequency to percentage scale
df_combined$Perc <- (df_combined$Freq / 1000) * 100

# Save the final combined data frame for use in figures
save(df_combined, file = "df_combined_Model_Based_Simulation_Sn.0.80_varying_b.Rdata")

# -------------------------------------------------------------------------------
# Figure 6
# Community Risk by Bag Size - Main Manuscript
# This section generates plots for community-level risk outcomes,
# based on model-based simulations under varying bag sizes (`b` values),
# using a diagnostic sensitivity (Sn) of 80%.
# -------------------------------------------------------------------------------

library(ggplot2)
library(ggpubr)
library(dplyr)

load("df_combined_Model_Based_Simulation_Sn.0.80_varying_b.Rdata")


# Figure 6 for Community Risk - Main Manuscript --------------#

df_nonzero <- df_combined[df_combined$Counts == "nonzero" , ]
df_zero <- df_combined[df_combined$Counts == "zero" , ]

# Keep only TP and FN
df_nonzero_TP_FN <- df_nonzero %>%
  filter(Outcome %in% c("TP","FN"))

# Step 1: compute total TP+FN per Tx and b
df_sums <- df_nonzero_TP_FN %>%
  group_by(Tx, b) %>%
  summarise(total_TP_FN = sum(Freq), .groups = "drop")

# Step 2: find max total TP+FN per Tx
df_max <- df_sums %>%
  group_by(Tx) %>%
  summarise(max_total = max(total_TP_FN, na.rm = TRUE), .groups = "drop")

# Step 3: join and rescale each Freq
df_fixed <- df_nonzero_TP_FN %>%
  left_join(df_sums, by = c("Tx", "b")) %>%
  left_join(df_max, by = "Tx") %>%
  mutate(Freq_fixed = ifelse(total_TP_FN > 0,
                             Freq / total_TP_FN * max_total,
                             0)) %>%
  select(b, Tx, Outcome, Freq_fixed)

# Check sums: should be identical across b for each Tx
df_fixed %>%
  group_by(Tx, b) %>%
  summarise(total = sum(Freq_fixed), .groups = "drop")

# Plot for nonzero counts
max_freq <- ceiling(max(df_fixed$Freq_fixed, na.rm = TRUE)) + 0.05

nonzero <- ggplot(df_fixed, aes(x = b, y = Freq_fixed, fill = Outcome)) +
  geom_col(position = "stack") +
  facet_wrap(~ Tx, nrow = 1) +
  scale_fill_manual(values = c("TP" = "blue", "TN" = "green", "FN" = "red")) +
  labs(x = "No. of sampled bags, b", title = paste("Tx > 0")) +
  scale_y_continuous(
    name = "Frequency (in 1000s)",
    breaks = seq(0, max_freq, by = 1),limits = c(0,max_freq)
  )+
  theme(
    strip.text = element_text(size = 8,face = "bold"),  # Reduce facet label font size to 8
    axis.text.x = element_text(size = 6, face = "bold")  # X-axis text
  )

# Plot for zero counts
max_freq_zero <- ceiling(max(df_zero$Freq, na.rm = TRUE)) + 0.05

df_zero$Freq <- round(df_zero$Freq)
zero <- ggplot(df_zero, aes(x = b, y = Freq, fill = Outcome)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("TP" = "blue", "TN" = "green", "FN" = "red")) +
  facet_wrap(~ Tx, nrow = 1) +
  labs(x = "b", title = paste("Tx = 0")) +
  scale_y_continuous(
    name = "Frequency (in 1000s)",
    breaks = seq(0, 1000, by = 100),limits = c(0,1000)
  )+
  theme(
    strip.text = element_text(size = 8,face = "bold"),  # Reduce facet label font size to 8
    axis.text.x = element_text(size = 6, face = "bold")  # X-axis text
  )

# Combine the zero and nonzero plots

    plot_m_b_13_20_26 <- annotate_figure(
      ggarrange(zero, nonzero, nrow = 1, common.legend = TRUE, legend = "top",
                widths = c(0.15, 0.85)))


    plot_m_b_13_20_26

ggsave("Test_Outcome_m5_b13_20_26_Imperfect_Manuscript.png",width = 10,height = 6,units = "in")

# --------------------------------------------------------------------------------
# Supplementary Table S5
# Community Risk assessment
# --------------------------------------------------------------------------------


# --------------------------------------------------------------------------------
# Leakage parameters: By Sampled bag size - b=13,20,26: k=0
# --------------------------------------------------------------------------------

load("submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect_Sn_80 <- NULL
load("submission\\Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect_Sn_80 <- NULL
load("submission\\Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect_Sn_80 <- NULL

LKG_TrBB_MH_00_m_5_b_13_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.7`)
LKG_TrBB_MH_00_m_5_b_13_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.8`)
LKG_TrBB_MH_00_m_5_b_13_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.9`)
LKG_TrBB_MH_00_m_5_b_13_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`0.95`)
LKG_TrBB_MH_00_m_5_b_13_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_00_Imperfect$`1`)

LKG_TrBB_MH_00_m_5_b_20_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect$`0.7`)
LKG_TrBB_MH_00_m_5_b_20_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect$`0.8`)
LKG_TrBB_MH_00_m_5_b_20_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect$`0.9`)
LKG_TrBB_MH_00_m_5_b_20_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect$`0.95`)
LKG_TrBB_MH_00_m_5_b_20_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_00_Imperfect$`1`)

LKG_TrBB_MH_00_m_5_b_26_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect$`0.7`)
LKG_TrBB_MH_00_m_5_b_26_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect$`0.8`)
LKG_TrBB_MH_00_m_5_b_26_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect$`0.9`)
LKG_TrBB_MH_00_m_5_b_26_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect$`0.95`)
LKG_TrBB_MH_00_m_5_b_26_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_00_Imperfect$`1`)

mean_LKG_df_MH_00_m_5_b_13 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_00_m_5_b_13_sn_0.7$E.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.8$E.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.9$E.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.95$E.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_00_m_5_b_13_sn_0.7$P.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.8$P.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.9$P.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.95$P.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_00_m_5_b_13_sn_0.7$Median.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.8$Median.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.9$Median.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.95$Median.Li,
    LKG_TrBB_MH_00_m_5_b_13_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_00_m_5_b_13_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_13_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_13_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_00_m_5_b_13_sn_0.7$P.Li)))
)

mean_LKG_df_MH_00_m_5_b_20 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_00_m_5_b_20_sn_0.7$E.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.8$E.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.9$E.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.95$E.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_00_m_5_b_20_sn_0.7$P.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.8$P.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.9$P.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.95$P.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_00_m_5_b_20_sn_0.7$Median.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.8$Median.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.9$Median.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.95$Median.Li,
    LKG_TrBB_MH_00_m_5_b_20_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_00_m_5_b_20_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_20_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_20_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_00_m_5_b_20_sn_0.7$P.Li)))
)

mean_LKG_df_MH_00_m_5_b_26 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_00_m_5_b_26_sn_0.7$E.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.8$E.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.9$E.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.95$E.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_00_m_5_b_26_sn_0.7$P.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.8$P.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.9$P.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.95$P.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_00_m_5_b_26_sn_0.7$Median.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.8$Median.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.9$Median.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.95$Median.Li,
    LKG_TrBB_MH_00_m_5_b_26_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_00_m_5_b_26_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_26_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_00_m_5_b_26_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_00_m_5_b_20_sn_0.7$P.Li)))
)


P_Li_MH_00_m_5_b_13 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_13, FUN = mean)
CI_P_Li_MH_00_m_5_b_13 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_13, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_00_m_5_b_13 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_13, FUN = mean)
CI_E_Li_MH_00_m_5_b_13 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_13, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))


P_Li_MH_00_m_5_b_20 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_20, FUN = mean)
CI_P_Li_MH_00_m_5_b_20 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_20, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_00_m_5_b_20 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_20, FUN = mean)
CI_E_Li_MH_00_m_5_b_20 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_20, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P_Li_MH_00_m_5_b_26 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_26, FUN = mean)
CI_P_Li_MH_00_m_5_b_26 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_26, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_00_m_5_b_26 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_26, FUN = mean)
CI_E_Li_MH_00_m_5_b_26 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_00_m_5_b_26, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P_Li_13_20_26_MH_00 <- cbind(P_Li_MH_00_m_5_b_13,CI_P_Li_MH_00_m_5_b_13$P.Li,P_Li_MH_00_m_5_b_20$P.Li,CI_P_Li_MH_00_m_5_b_20$P.Li,P_Li_MH_00_m_5_b_26$P.Li,CI_P_Li_MH_00_m_5_b_26$P.Li)
P_Li_13_20_26_MH_00[, sapply(P_Li_13_20_26_MH_00, is.numeric)] <-
  round(P_Li_13_20_26_MH_00[, sapply(P_Li_13_20_26_MH_00, is.numeric)] * 100,2)
names(P_Li_13_20_26_MH_00) <- c("Sensitivity","b_13","b_13_ll","b_13_ul","b_20","b_20_ll","b_20_ul","b_26","b_26_ll","b_26_ul")
E_Li_13_20_26_MH_00 <- cbind(E_Li_MH_00_m_5_b_13,CI_E_Li_MH_00_m_5_b_13$E.Li,E_Li_MH_00_m_5_b_20$E.Li,CI_E_Li_MH_00_m_5_b_20$E.Li,E_Li_MH_00_m_5_b_26$E.Li,CI_E_Li_MH_00_m_5_b_26$E.Li)
E_Li_13_20_26_MH_00[, sapply(E_Li_13_20_26_MH_00, is.numeric)] <-
  round(E_Li_13_20_26_MH_00[, sapply(E_Li_13_20_26_MH_00, is.numeric)],2)
names(E_Li_13_20_26_MH_00) <- c("Sensitivity","b_13","b_13_ll","b_13_ul","b_20","b_20_ll","b_20_ul","b_26","b_26_ll","b_26_ul")


P_E_Li_13_20_26_MH_00 <- rbind(P_Li_13_20_26_MH_00,E_Li_13_20_26_MH_00)

write.csv(P_E_Li_13_20_26_MH_00,file="P_E_Li_13_20_26_MH_00.csv")

P_E_Li_13_20_26_MH_00 <- read.csv("P_E_Li_13_20_26_MH_00.csv",header = TRUE)

print(xtable(P_E_Li_13_20_26_MH_00),rownames=FALSE)

# --------------------------------------------------------------------------------
# Leakage parameters: By Sampled bag size - b=13,20,26: k=0.01
# --------------------------------------------------------------------------------


load("submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect_Sn_80 <- NULL
load("submission\\Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect_Sn_80 <- NULL
load("submission\\Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80.Rdata")
Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect_Sn_80 <- NULL

LKG_TrBB_MH_01_m_5_b_13_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.7`)
LKG_TrBB_MH_01_m_5_b_13_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.8`)
LKG_TrBB_MH_01_m_5_b_13_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.9`)
LKG_TrBB_MH_01_m_5_b_13_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`0.95`)
LKG_TrBB_MH_01_m_5_b_13_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_01_Imperfect$`1`)

LKG_TrBB_MH_01_m_5_b_20_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect$`0.7`)
LKG_TrBB_MH_01_m_5_b_20_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect$`0.8`)
LKG_TrBB_MH_01_m_5_b_20_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect$`0.9`)
LKG_TrBB_MH_01_m_5_b_20_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect$`0.95`)
LKG_TrBB_MH_01_m_5_b_20_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_01_Imperfect$`1`)

LKG_TrBB_MH_01_m_5_b_26_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect$`0.7`)
LKG_TrBB_MH_01_m_5_b_26_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect$`0.8`)
LKG_TrBB_MH_01_m_5_b_26_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect$`0.9`)
LKG_TrBB_MH_01_m_5_b_26_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect$`0.95`)
LKG_TrBB_MH_01_m_5_b_26_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_01_Imperfect$`1`)


mean_LKG_df_MH_01_m_5_b_13 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_01_m_5_b_13_sn_0.7$E.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.8$E.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.9$E.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.95$E.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_01_m_5_b_13_sn_0.7$P.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.8$P.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.9$P.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.95$P.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_01_m_5_b_13_sn_0.7$Median.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.8$Median.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.9$Median.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.95$Median.Li,
    LKG_TrBB_MH_01_m_5_b_13_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_01_m_5_b_13_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_13_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_13_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_01_m_5_b_13_sn_0.7$P.Li)))
)

mean_LKG_df_MH_01_m_5_b_20 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_01_m_5_b_20_sn_0.7$E.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.8$E.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.9$E.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.95$E.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_01_m_5_b_20_sn_0.7$P.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.8$P.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.9$P.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.95$P.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_01_m_5_b_20_sn_0.7$Median.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.8$Median.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.9$Median.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.95$Median.Li,
    LKG_TrBB_MH_01_m_5_b_20_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_01_m_5_b_20_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_20_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_20_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_01_m_5_b_20_sn_0.7$P.Li)))
)


mean_LKG_df_MH_01_m_5_b_26 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_01_m_5_b_26_sn_0.7$E.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.8$E.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.9$E.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.95$E.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_01_m_5_b_26_sn_0.7$P.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.8$P.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.9$P.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.95$P.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_01_m_5_b_26_sn_0.7$Median.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.8$Median.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.9$Median.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.95$Median.Li,
    LKG_TrBB_MH_01_m_5_b_26_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_01_m_5_b_26_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_26_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_01_m_5_b_26_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_01_m_5_b_20_sn_0.7$P.Li)))
)

P_Li_MH_01_m_5_b_13 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_13, FUN = mean)
CI_P_Li_MH_01_m_5_b_13 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_13, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_01_m_5_b_13 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_13, FUN = mean)
CI_E_Li_MH_01_m_5_b_13 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_13, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))


P_Li_MH_01_m_5_b_20 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_20, FUN = mean)
CI_P_Li_MH_01_m_5_b_20 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_20, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_01_m_5_b_20 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_20, FUN = mean)
CI_E_Li_MH_01_m_5_b_20 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_20, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P_Li_MH_01_m_5_b_26 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_26, FUN = mean)
CI_P_Li_MH_01_m_5_b_26 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_26, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_01_m_5_b_26 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_26, FUN = mean)
CI_E_Li_MH_01_m_5_b_26 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_01_m_5_b_26, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P_Li_13_20_26_MH_01 <- cbind(P_Li_MH_01_m_5_b_13,CI_P_Li_MH_01_m_5_b_13$P.Li,P_Li_MH_01_m_5_b_20$P.Li,CI_P_Li_MH_01_m_5_b_20$P.Li,P_Li_MH_01_m_5_b_26$P.Li,CI_P_Li_MH_01_m_5_b_26$P.Li)
P_Li_13_20_26_MH_01[, sapply(P_Li_13_20_26_MH_01, is.numeric)] <-
  round(P_Li_13_20_26_MH_01[, sapply(P_Li_13_20_26_MH_01, is.numeric)] * 100,2)
names(P_Li_13_20_26_MH_01) <- c("Sensitivity","b_13","b_13_ll","b_13_ul","b_20","b_20_ll","b_20_ul","b_26","b_26_ll","b_26_ul")
E_Li_13_20_26_MH_01 <- cbind(E_Li_MH_01_m_5_b_13,CI_E_Li_MH_01_m_5_b_13$E.Li,E_Li_MH_01_m_5_b_20$E.Li,CI_E_Li_MH_01_m_5_b_20$E.Li,E_Li_MH_01_m_5_b_26$E.Li,CI_E_Li_MH_01_m_5_b_26$E.Li)
E_Li_13_20_26_MH_01[, sapply(E_Li_13_20_26_MH_01, is.numeric)] <-
  round(E_Li_13_20_26_MH_01[, sapply(E_Li_13_20_26_MH_01, is.numeric)],2)
names(E_Li_13_20_26_MH_01) <- c("Sensitivity","b_13","b_13_ll","b_13_ul","b_20","b_20_ll","b_20_ul","b_26","b_26_ll","b_26_ul")


P_E_Li_13_20_26_MH_01 <- rbind(P_Li_13_20_26_MH_01,E_Li_13_20_26_MH_01)


write.csv(P_E_Li_13_20_26_MH_01,file="P_E_Li_13_20_26_MH_01.csv")


P_E_Li_13_20_26_MH_01 <- read.csv("P_E_Li_13_20_26_MH_01.csv",header=TRUE)

library(xtable)
print(xtable(P_E_Li_13_20_26_MH_01),rownames=FALSE)

# --------------------------------------------------------------------------------
# Leakage parameters: By Sampled bag size - b=13,20,26: k=0.02
# --------------------------------------------------------------------------------

load("submission\\Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80.Rdata")
load("submission\\Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect_Sn_80.Rdata")
load("submission\\Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect_Sn_80.Rdata")

Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect_Sn_80 <- NULL
Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect_Sn_80 <- NULL
Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect <- Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect_Sn_80
Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect_Sn_80 <- NULL

LKG_TrBB_MH_02_m_5_b_13_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.7`)
LKG_TrBB_MH_02_m_5_b_13_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.8`)
LKG_TrBB_MH_02_m_5_b_13_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.9`)
LKG_TrBB_MH_02_m_5_b_13_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`0.95`)
LKG_TrBB_MH_02_m_5_b_13_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b13_MH_02_Imperfect$`1`)

LKG_TrBB_MH_02_m_5_b_20_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect$`0.7`)
LKG_TrBB_MH_02_m_5_b_20_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect$`0.8`)
LKG_TrBB_MH_02_m_5_b_20_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect$`0.9`)
LKG_TrBB_MH_02_m_5_b_20_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect$`0.95`)
LKG_TrBB_MH_02_m_5_b_20_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b20_MH_02_Imperfect$`1`)

LKG_TrBB_MH_02_m_5_b_26_sn_0.7 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect$`0.7`)
LKG_TrBB_MH_02_m_5_b_26_sn_0.8 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect$`0.8`)
LKG_TrBB_MH_02_m_5_b_26_sn_0.9 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect$`0.9`)
LKG_TrBB_MH_02_m_5_b_26_sn_0.95 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect$`0.95`)
LKG_TrBB_MH_02_m_5_b_26_sn_1 <- compute_simulated_leakage(Sim_Tx_TrBB_Prawn_m5_b26_MH_02_Imperfect$`1`)


# Combine leakage values in dataframe

mean_LKG_df_MH_02_m_5_b_13 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_02_m_5_b_13_sn_0.7$E.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.8$E.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.9$E.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.95$E.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_02_m_5_b_13_sn_0.7$P.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.8$P.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.9$P.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.95$P.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_02_m_5_b_13_sn_0.7$Median.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.8$Median.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.9$Median.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.95$Median.Li,
    LKG_TrBB_MH_02_m_5_b_13_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_02_m_5_b_13_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_13_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_13_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_02_m_5_b_13_sn_0.7$P.Li)))
)

mean_LKG_df_MH_02_m_5_b_20 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_02_m_5_b_20_sn_0.7$E.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.8$E.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.9$E.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.95$E.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_02_m_5_b_20_sn_0.7$P.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.8$P.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.9$P.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.95$P.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_02_m_5_b_20_sn_0.7$Median.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.8$Median.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.9$Median.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.95$Median.Li,
    LKG_TrBB_MH_02_m_5_b_20_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_02_m_5_b_20_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_20_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_20_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_02_m_5_b_20_sn_0.7$P.Li)))
)



mean_LKG_df_MH_02_m_5_b_26 <- data.frame(
  E.Li = c(
    LKG_TrBB_MH_02_m_5_b_26_sn_0.7$E.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.8$E.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.9$E.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.95$E.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_1$E.Li
  ),
  P.Li = c(
    LKG_TrBB_MH_02_m_5_b_26_sn_0.7$P.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.8$P.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.9$P.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.95$P.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_1$P.Li
  ),
  Median.Li = c(
    LKG_TrBB_MH_02_m_5_b_26_sn_0.7$Median.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.8$Median.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.9$Median.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.95$Median.Li,
    LKG_TrBB_MH_02_m_5_b_26_sn_1$Median.Li
  ),
  Pseudo.E.Li.median = c(
    LKG_TrBB_MH_02_m_5_b_26_sn_0.7$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.8$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.9$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_26_sn_0.95$Pseudo.E.Li.median,
    LKG_TrBB_MH_02_m_5_b_26_sn_1$Pseudo.E.Li.median
  ),
  sensitivity = factor(rep(c(0.7, 0.8, 0.9, 0.95, 1), each = length(LKG_TrBB_MH_02_m_5_b_20_sn_0.7$P.Li)))
)


# --------------------------------------------------------------------------------
# Calculate mean and 95% Credible interval by sensitivity, k and bag size
# --------------------------------------------------------------------------------



P_Li_MH_02_m_5_b_13 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_13, FUN = mean)
CI_P_Li_MH_02_m_5_b_13 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_13, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_02_m_5_b_13 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_13, FUN = mean)
CI_E_Li_MH_02_m_5_b_13 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_13, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))


P_Li_MH_02_m_5_b_20 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_20, FUN = mean)
CI_P_Li_MH_02_m_5_b_20 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_20, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_02_m_5_b_20 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_20, FUN = mean)
CI_E_Li_MH_02_m_5_b_20 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_20, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P_Li_MH_02_m_5_b_26 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_26, FUN = mean)
CI_P_Li_MH_02_m_5_b_26 <- aggregate(P.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_26, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))
E_Li_MH_02_m_5_b_26 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_26, FUN = mean)
CI_E_Li_MH_02_m_5_b_26 <- aggregate(E.Li ~ sensitivity, data = mean_LKG_df_MH_02_m_5_b_26, FUN = function(x) quantile(x, probs = c(0.025, 0.975)))

P_Li_13_20_26_MH_02 <- cbind(P_Li_MH_02_m_5_b_13,CI_P_Li_MH_02_m_5_b_13$P.Li,P_Li_MH_02_m_5_b_20$P.Li,CI_P_Li_MH_02_m_5_b_20$P.Li,P_Li_MH_02_m_5_b_26$P.Li,CI_P_Li_MH_02_m_5_b_26$P.Li)
P_Li_13_20_26_MH_02[, sapply(P_Li_13_20_26_MH_02, is.numeric)] <-
  round(P_Li_13_20_26_MH_02[, sapply(P_Li_13_20_26_MH_02, is.numeric)] * 102,2)
names(P_Li_13_20_26_MH_02) <- c("Sensitivity","b_13","b_13_ll","b_13_ul","b_20","b_20_ll","b_20_ul","b_26","b_26_ll","b_26_ul")
E_Li_13_20_26_MH_02 <- cbind(E_Li_MH_02_m_5_b_13,CI_E_Li_MH_02_m_5_b_13$E.Li,E_Li_MH_02_m_5_b_20$E.Li,CI_E_Li_MH_02_m_5_b_20$E.Li,E_Li_MH_02_m_5_b_26$E.Li,CI_E_Li_MH_02_m_5_b_26$E.Li)
E_Li_13_20_26_MH_02[, sapply(E_Li_13_20_26_MH_02, is.numeric)] <-
  round(E_Li_13_20_26_MH_02[, sapply(E_Li_13_20_26_MH_02, is.numeric)],2)
names(E_Li_13_20_26_MH_02) <- c("Sensitivity","b_13","b_13_ll","b_13_ul","b_20","b_20_ll","b_20_ul","b_26","b_26_ll","b_26_ul")


P_E_Li_13_20_26_MH_02 <- rbind(P_Li_13_20_26_MH_02,E_Li_13_20_26_MH_02)

write.csv(P_E_Li_13_20_26_MH_02,file="P_E_Li_13_20_26_MH_02.csv")

P_E_Li_13_20_26_MH_02 <- read.csv("P_E_Li_13_20_26_MH_02.csv",header = TRUE)

library(xtable)
print(xtable(P_E_Li_13_20_26_MH_02),rownames=FALSE)

# --------------------------------------------------------------------------------
# Final Results for Supplementary Table S5
# --------------------------------------------------------------------------------

P_E_Li_13_20_26_MH_00 <- read.csv("P_E_Li_13_20_26_MH_00.csv",header = TRUE)
P_E_Li_13_20_26_MH_01 <- read.csv("P_E_Li_13_20_26_MH_01.csv",header = TRUE)
P_E_Li_13_20_26_MH_02 <- read.csv("P_E_Li_13_20_26_MH_02.csv",header = TRUE)

P_E_Li_13_20_26_MH_00_01_02 <- rbind(P_E_Li_13_20_26_MH_00,
                                     P_E_Li_13_20_26_MH_01,
                                     P_E_Li_13_20_26_MH_02)
library(xtable)
print(xtable(P_E_Li_13_20_26_MH_00_01_02),rownames=FALSE)
