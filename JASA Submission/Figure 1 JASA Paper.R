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

#-------------------------------------------------------------@
# Figure 1: JASA PAper
# Supplementary Figure 2: JASA PAper
# Data: Deidentified Frozen Seafood Data
# Exact BB Model: MLE method and PLLF
# Using Generailzed Function
# Target Parameters: PLLF of Sensitivity and Specificity
#-------------------------------------------------------------@

#-------------------------------------------------------------
# Deidentified Frozen Seafood Data
#-------------------------------------------------------------@

df.seafood <- read.csv("submission//deidentified_frozen_seafood.csv")
df.seafood$Yfac <- factor(df.seafood$numPos,levels=c(0:13))
x <- as.numeric(names(table(df.seafood$Yfac)))
freq <- as.numeric(table(df.seafood$Yfac))
size <- 13 # Number of bags (b)
seafood.data <- data.frame(ty=df.seafood$Yfac,n=rep(size,length(df.seafood$Yfac)))
seafood.data$ID <- c(1:dim(seafood.data)[1])
seafood.data$ty <- as.numeric(paste(seafood.data$ty))
seafood.data <- seafood.data[order(seafood.data$ty),]

summ.seafood.data <- as.data.frame(table(seafood.data$ty))
colnames(summ.seafood.data) <- c("ty","freq")
summ.seafood.data$ty <- as.numeric(paste(summ.seafood.data$ty))


#-------------------------------------------------------------
# PLLF: Specificity under imperfect sensitivity (0.8)
#-------------------------------------------------------------

# Function to run the optimization for a given test specificity value
store_specificity_values_prawn_imperfect <- function(specificity) {

  # Run the optimization using the log-likelihood function 'loglik_group_trbb'
  optim_result <- optim(
    par = c(0, 0),                     # Initial parameter guesses
    fn = loglik_group_trbb,            # Log-likelihood function
    ty = summ.seafood.data$ty,         # Observed data: positive tests per group
    freq = summ.seafood.data$freq,     # Frequency of each observation
    b = 13,                            # Number of sampled groups per consignment
    m = 5,                             # Number of items per group
    M = 40,                            # Total number of items per consignment
    theta = Inf,                       # Overdispersion parameter to consider No clustering
    R = 1e4,                           # Number of static quantile from Beta distribution
    sensitivity = 0.8,                 # Fixed test sensitivity
    specificity = specificity,         # Input test specificity
    cutoff = 0,                        # Detection cutoff
    deviance = FALSE,                  # Whether to return deviance
    control = list(reltol = 1e-12, fnscale = -1),  # Optim control: maximize LL
    hessian = FALSE                    # No Hessian matrix needed
  )

  # Extract and compute model parameters
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$cutoff <- 0
  optim_result$mu.TrBB <- optim_result$alpha / (optim_result$alpha + optim_result$beta)
  optim_result$sensitivity <- 0.8
  optim_result$specificity <- specificity

  # Return results as a data frame
  return(data.frame(
    cutoff = optim_result$cutoff,
    alpha = optim_result$alpha,
    beta = optim_result$beta,
    mu.BB = optim_result$mu.TrBB,
    sensitivity = optim_result$sensitivity,
    specificity = optim_result$specificity,
    value = optim_result$value,
    convergence = optim_result$convergence
  ))
}

# Define a sequence of specificity values to test
specificity <- seq(0.9989, 1, by = 0.000001)

# Apply the function to each specificity value using lapply()
# Apply the function with error handling
PLLF_Specificity_imperfect_Sensitivity_Prawn <-
  lapply(specificity, store_specificity_values_prawn_imperfect)

# Apply the function with error handling
PLLF_Specificity_imperfect_Sensitivity_Prawn <- lapply(specificity, function(spcs) {
  tryCatch({
    # Run the optimization for this specificity value
    store_specificity_values_prawn_imperfect(spcs)
  },
  error = function(e) {
    # Handle any error without stopping the loop
    message("Error for specificity = ", spcs, ": ", e$message)
    return(NA)  # or NULL, depending on your downstream needs
  })
})

# Remove any failed (NULL) results from the list
PLLF_Specificity_imperfect_Sensitivity_Prawn <-
  Filter(Negate(is.null), PLLF_Specificity_imperfect_Sensitivity_Prawn)

# Combine the list of successful results into one data frame
PLLF_Specificity_imperfect_Sensitivity_Prawn <-
  do.call(rbind, PLLF_Specificity_imperfect_Sensitivity_Prawn)

# Save the combined results to an .Rdata file
save(
  PLLF_Specificity_imperfect_Sensitivity_Prawn,
  file = "PLLF_Specificity_imperfect_Sensitivity_Prawn_Re4.Rdata"
)



#-------------------------------------------------------------
# PLLF: Sensitivity  under perfect Specificity (1.0)
#-------------------------------------------------------------

# Function to run the optimization for a given test sensitivity value
# (Assumes perfect specificity = 1.0)
store_sensitivity_values_prawn_perfect <- function(sensitivity) {

  # Run the optimization using the log-likelihood function 'loglik_group_trbb'
  optim_result <- optim(
    par = c(0, 0),                     # Initial parameter guesses
    fn = loglik_group_trbb,            # Log-likelihood function
    ty = summ.seafood.data$ty,         # Observed data: positive tests per group
    freq = summ.seafood.data$freq,     # Frequency of each observation
    b = 13,                            # Number of sampled groups per consignment
    m = 5,                             # Number of items per group
    M = 40,                            # Total number of items per consignment
    theta = Inf,                       # Overdispersion parameter to consider No clustering
    R = 1e4,                           # Number of static quantile from Beta distribution
    sensitivity = sensitivity,         # Input sensitivity value
    specificity = 1.0,                 # Fixed perfect specificity
    cutoff = 0,                        # Detection cutoff
    deviance = FALSE,                  # Whether to return deviance
    control = list(reltol = 1e-12, fnscale = -1),  # Optimization control
    hessian = FALSE                    # No Hessian matrix needed
  )

  # Extract and compute model parameters
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$cutoff <- 0
  optim_result$mu.TrBB <- optim_result$alpha / (optim_result$alpha + optim_result$beta)
  optim_result$sensitivity <- sensitivity
  optim_result$specificity <- 1.0

  # Return results as a data frame
  return(data.frame(
    cutoff = optim_result$cutoff,
    alpha = optim_result$alpha,
    beta = optim_result$beta,
    mu.BB = optim_result$mu.TrBB,
    sensitivity = optim_result$sensitivity,
    specificity = optim_result$specificity,
    value = optim_result$value,
    convergence = optim_result$convergence
  ))
}

# Define a sequence of sensitivity values to test
sensitivity <- seq(0.45, 1, by = 0.0001)

# Apply the function to each sensitivity value using lapply()
PLLF_Sensitivity_perfect_Specificity_Prawn <-
  lapply(sensitivity, store_sensitivity_values_prawn_perfect)
# Apply the function with error handling
PLLF_Sensitivity_perfect_Specificity_Prawn <- lapply(sensitivity, function(snsv) {
  tryCatch({
    # Run the optimization for this specificity value
    store_sensitivity_values_prawn_perfect(snsv)
  },
  error = function(e) {
    # Handle any error without stopping the loop
    message("Error for sensitivity = ", snsv, ": ", e$message)
    return(NA)  # or NULL, depending on your downstream needs
  })
})

# Combine all results into a single data frame
PLLF_Sensitivity_perfect_Specificity_Prawn <-
  do.call(rbind, PLLF_Sensitivity_perfect_Specificity_Prawn)

# Save the combined results to an .Rdata file
save(
  PLLF_Sensitivity_perfect_Specificity_Prawn,
  file = "PLLF_Sensitivity_perfect_Specificity_Prawn_Re4.Rdata"
)

#-------------------------------------------------------------
# PLLF: Sensitivity and specificity
# Figure 1: Manuscript
#-------------------------------------------------------------

load("PLLF_Specificity_imperfect_Sensitivity_Prawn_Re4.Rdata") # Using Generailzed function
load("PLLF_Sensitivity_perfect_Specificity_Prawn_Re4.Rdata")   # Using Generailzed function

df.pllf.sn.sp.1.prawn.0.45<-PLLF_Sensitivity_perfect_Specificity_Prawn
df.pllf.sn.sp.1.prawn.0.45$logL <- df.pllf.sn.sp.1.prawn.0.45$value
df.pllf.sn.sp.1.prawn.0.45$dif.logL.ll <- abs(df.pllf.sn.sp.1.prawn.0.45$logL - (max(df.pllf.sn.sp.1.prawn.0.45$logL,na.rm=T) -3.84/2))
ll.sensitivity <- df.pllf.sn.sp.1.prawn.0.45$sensitivity[which(df.pllf.sn.sp.1.prawn.0.45$dif.logL.ll==min(df.pllf.sn.sp.1.prawn.0.45$dif.logL.ll))]
mle.sensitivity <- df.pllf.sn.sp.1.prawn.0.45$sensitivity[which(df.pllf.sn.sp.1.prawn.0.45$logL==max(df.pllf.sn.sp.1.prawn.0.45$logL,na.rm=T))]
ul.sensitivity <- 1

pllf.CI.estimate.sensitivity.95 <- pllf.CI.estimate(df.pllf.sn.sp.1.prawn.0.45$sensitivity,df.pllf.sn.sp.1.prawn.0.45$logL,level = 0.95)

df.pllf.sp.sn.1.prawn.0.999<-PLLF_Specificity_imperfect_Sensitivity_Prawn
df.pllf.sp.sn.1.prawn.0.999$logL <- df.pllf.sp.sn.1.prawn.0.999$value
df.pllf.sp.sn.1.prawn.0.999$dif.logL.ll <- abs(df.pllf.sp.sn.1.prawn.0.999$logL - (max(df.pllf.sp.sn.1.prawn.0.999$logL,na.rm=T) -3.84/2))
ll.specificity <- df.pllf.sp.sn.1.prawn.0.999$specificity[which(df.pllf.sp.sn.1.prawn.0.999$dif.logL.ll==min(df.pllf.sp.sn.1.prawn.0.999$dif.logL.ll))]
mle.specificity <- df.pllf.sp.sn.1.prawn.0.999$specificity[which(df.pllf.sp.sn.1.prawn.0.999$logL==max(df.pllf.sp.sn.1.prawn.0.999$logL,na.rm=T))]
ul.specificity <- 1

pllf.CI.estimate.specificity.95 <- pllf.CI.estimate(df.pllf.sp.sn.1.prawn.0.999$specificity,df.pllf.sp.sn.1.prawn.0.999$logL,level = 0.95)


# ---------------------------------------------#
# Figure 1 : On Original scael of logL
# ---------------------------------------------#
png(file="pllf CI sensitivity specifity.png", width = 8,height = 4,units = "in",res = 300)

par(mfrow=c(1,2), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 4, 4, 2))  # Adjusted margins

range_ll <- range(df.pllf.sn.sp.1.prawn.0.45$logL)

plot(df.pllf.sn.sp.1.prawn.0.45$sensitivity[],df.pllf.sn.sp.1.prawn.0.45$logL,typ="l",xlab="Sensitivity",ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1) # ,main="logL under perfect specificity"
axis(side = 2, at = seq(round(range_ll[1],1), round(range_ll[2],1), by = 0.10),las = 1, asp = 1)
axis(side = 1, at = seq(min(df.pllf.sn.sp.1.prawn.0.45$sensitivity), max(df.pllf.sn.sp.1.prawn.0.45$sensitivity), by = 0.05),las = 2)
abline(h=max(df.pllf.sn.sp.1.prawn.0.45$logL,na.rm=T),col=2,lty=3)
abline(h=max(df.pllf.sn.sp.1.prawn.0.45$logL,na.rm=T)-3.84/2,col=3)
abline(v=pllf.CI.estimate.sensitivity.95$mle,col=4,lty=1)
abline(v=pllf.CI.estimate.sensitivity.95$mle.ll,col=4,lty=2)
abline(v=pllf.CI.estimate.sensitivity.95$mle.ul,col=4,lty=2)
# legend("bottomright",legend=c(paste("MLE=",mle.sensitivity),paste("LL=",ll.sensitivity),paste("UL=",ul.sensitivity)),bty = "n")


# Identify the x-value corresponding to the maximum y-value
max_x <- pllf.CI.estimate.sensitivity.95$mle
max_y <- max(df.pllf.sn.sp.1.prawn.0.45$logL,na.rm=T)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
text(x=max_x-0.01, y=(max_y-1.5), labels=round(max_x, 5), pos=3, col=2, cex=0.9, srt = 90)
text(x=max_x+0.25, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.95$mle.ul, 5), pos=3, col=6, cex=0.9, srt = 90)
text(x=max_x-0.27, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.95$mle.ll, 5), pos=3, col=6, cex=0.9, srt = 90)



range_ll <- range(df.pllf.sp.sn.1.prawn.0.999$logL)
plot(df.pllf.sp.sn.1.prawn.0.999$specificity,df.pllf.sp.sn.1.prawn.0.999$logL,typ="l",xlab="Specificity",ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1) # ,main="logL under perfect sensitivity"
axis(side = 2, at = seq(round(range_ll[1],1), round(range_ll[2]+1,1), by = 2),las = 1, asp = 1)
axis(side = 1, at = seq(min(df.pllf.sn.sp.1.prawn.0.45$sensitivity), max(df.pllf.sn.sp.1.prawn.0.45$sensitivity), by = 0.0001),las = 2)
abline(h=max(df.pllf.sp.sn.1.prawn.0.999$logL),col=2,lty=1)
abline(h=max(df.pllf.sp.sn.1.prawn.0.999$logL)-3.84/2,col=3,lty=1)
abline(v=mle.specificity,col=4,lty=1)
abline(v=ll.specificity,col=4,lty=2)
# legend("bottomright",legend=c(paste("MLE=",mle.specificity),paste("LL=",ll.specificity),paste("UL=",ul.specificity)),bty = "n")
# legend("bottomright",legend=c("MLE=1.0","LL=0.999881"),bty = "n")

# Identify the x-value corresponding to the maximum y-value
max_x <- pllf.CI.estimate.specificity.95$mle
max_y <- max(df.pllf.sp.sn.1.prawn.0.999$logL,na.rm=T)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
text(x=max_x-0.00002, y=(max_y-15), labels=round(max_x, 5), pos=3, col=2, cex=0.9, srt = 90)
# text(x=max_x+0.25, y=(max_y-2), labels=round(pllf.CI.estimate.specificity.95$mle.ul, 5), pos=3, col=6, cex=0.9, srt = 90)
text(x=max_x-0.00015, y=(max_y-15), labels=round(pllf.CI.estimate.specificity.95$mle.ll, 5), pos=3, col=6, cex=0.9, srt = 90)


dev.off()

# ---------------------------------------------#
# Figure 1 : On Re-scaled logL  (logL - maxlogL)
# ---------------------------------------------#

png(file="pllf CI sensitivity specifity rescaled.png", width = 8,height = 4,units = "in",res = 300)

par(mfrow=c(1,2), mgp=c(3.25, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 5, 4, 2))  # Adjusted margins

df.pllf.sn.sp.1.prawn.0.45$logL_rescaled <- df.pllf.sn.sp.1.prawn.0.45$logL - max(df.pllf.sn.sp.1.prawn.0.45$logL,na.rm=T)

range_ll <- range(df.pllf.sn.sp.1.prawn.0.45$logL_rescaled)

step <- -0.25
ymax <- 0
ymin <- floor(range_ll[1]*10)/10  # round down to nearest 0.1
y_ticks <- seq(from = ymax, to = ymin, by = step)

plot(df.pllf.sn.sp.1.prawn.0.45$sensitivity,
     df.pllf.sn.sp.1.prawn.0.45$logL_rescaled,
     typ="l",xlab="Sensitivity",
     ylab="logL - max(logL)",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1,
     ylim = c(range_ll[1], 0)) # ,main="logL under perfect specificity"
axis(side = 2, at = y_ticks,las = 1,tcl = -0.3)
axis(side = 1, at = seq(min(df.pllf.sn.sp.1.prawn.0.45$sensitivity), max(df.pllf.sn.sp.1.prawn.0.45$sensitivity), by = 0.05),
     las = 2,tcl = -0.3)
abline(h=max(df.pllf.sn.sp.1.prawn.0.45$logL_rescaled,na.rm=T),col=2,lty=2)
abline(h=max(df.pllf.sn.sp.1.prawn.0.45$logL_rescaled,na.rm=T)-3.84/2,col=4,lty=3)
abline(h=max(df.pllf.sn.sp.1.prawn.0.45$logL_rescaled,na.rm=T)-2.71/2,col=5,lty=6)
#abline(v=pllf.CI.estimate.sensitivity.95$mle,col=4,lty=1)
#abline(v=pllf.CI.estimate.sensitivity.95$mle.ll,col=4,lty=2)
#abline(v=pllf.CI.estimate.sensitivity.95$mle.ul,col=4,lty=2)
# legend("bottomright",legend=c(paste("MLE=",mle.sensitivity),paste("LL=",ll.sensitivity),paste("UL=",ul.sensitivity)),bty = "n")

legend("bottomright",legend=c("0.95","0.90"),col=c(4,5),lty=c(3,6),bty="n",lwd=c(2,2))


# Identify the x-value corresponding to the maximum y-value
#max_x <- pllf.CI.estimate.sensitivity.95$mle
#max_y <- max(df.pllf.sn.sp.1.prawn.0.45$logL_rescaled,na.rm=T)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
#text(x=max_x-0.01, y=(max_y-1.5), labels=round(max_x, 5), pos=3, col=2, cex=0.9, srt = 90)
#text(x=max_x+0.25, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.95$mle.ul, 5), pos=3, col=6, cex=0.9, srt = 90)
#text(x=max_x-0.27, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.95$mle.ll, 5), pos=3, col=6, cex=0.9, srt = 90)


df.pllf.sp.sn.1.prawn.0.999$logL_rescaled <- df.pllf.sp.sn.1.prawn.0.999$logL - max(df.pllf.sp.sn.1.prawn.0.999$logL,na.rm=T)
range_ll <- range(df.pllf.sp.sn.1.prawn.0.999$logL_rescaled)

step <- -2.5
ymax <- 0
ymin <- floor(range_ll[1]*10)/10  # round down to nearest 0.1
y_ticks <- seq(from = ymax, to = ymin, by = step)


plot(df.pllf.sp.sn.1.prawn.0.999$specificity,
     df.pllf.sp.sn.1.prawn.0.999$logL_rescaled,typ="l",
     xlab="Specificity",ylab="logL - max(logL)",
     yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1,
     ylim = c(range_ll[1], 0)) # ,main="logL under perfect sensitivity"
axis(side = 2, at = y_ticks,las = 1,tcl = -0.3,cex.axis=0.9)
axis(side = 1, at = seq(min(df.pllf.sp.sn.1.prawn.0.999$specificity), max(df.pllf.sp.sn.1.prawn.0.999$specificity), by = 0.0001),
     las = 2,tcl = -0.3,cex.axis=0.9)
abline(h=max(df.pllf.sp.sn.1.prawn.0.999$logL_rescaled),col=2,lty=2)
abline(h=max(df.pllf.sp.sn.1.prawn.0.999$logL_rescaled)-3.84/2,col=4,lty=3)
abline(h=max(df.pllf.sp.sn.1.prawn.0.999$logL_rescaled)-2.71/2,col=5,lty=4)
# abline(v=mle.specificity,col=4,lty=1)
# abline(v=ll.specificity,col=4,lty=2)
# legend("bottomright",legend=c(paste("MLE=",mle.specificity),paste("LL=",ll.specificity),paste("UL=",ul.specificity)),bty = "n")
# legend("bottomright",legend=c("MLE=1.0","LL=0.999881"),bty = "n")

# Identify the x-value corresponding to the maximum y-value
#max_x <- pllf.CI.estimate.specificity.95$mle
#max_y <- max(df.pllf.sp.sn.1.prawn.0.999$logL_rescaled,na.rm=T)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
#text(x=max_x-0.00002, y=(max_y-15), labels=round(max_x, 5), pos=3, col=2, cex=0.9, srt = 90)
# text(x=max_x+0.25, y=(max_y-2), labels=round(pllf.CI.estimate.specificity.95$mle.ul, 5), pos=3, col=6, cex=0.9, srt = 90)
#text(x=max_x-0.00015, y=(max_y-15), labels=round(pllf.CI.estimate.specificity.95$mle.ll, 5), pos=3, col=6, cex=0.9, srt = 90)

legend("bottomright",legend=c("0.95","0.90"),col=c(4,5),lty=c(3,6),bty="n",lwd=c(2,2))

dev.off()


# ---------------------------------------------------------
# Supplementary Material Analysis
# PLLF for Sensitivity: Without restriction on Specificity
# PLLF for Specificity: Without restriction on sensitivity
# Using Cumulative Beta Distribution: loglik_group_bb
# ---------------------------------------------------------


# ---------------------------------------------------------#
# PLLF for specificity: Without restriction on Sensitivity
# ---------------------------------------------------------#

# MLE Estimates: alpha, beta and sensitivity
o.theta.inf.prawn.sp.sn <- optim(c(0,0,0),loglik_group_bb,ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,
                                 theta=Inf,m=5,M=40,specificity=1,deviance=FALSE,control=list(factr=1e-12,fnscale=-1,trace=T),
                                 R=1e4,hessian=FALSE,method = "L-BFGS-B",lower = c(log(0.0001),log(0.0001),log(0.50)),upper = c(Inf,Inf,0))

# PLLF Calculation: Specificity without restricting Sensitivity
# Function to run the optimization for a given specificity value
store_specificity_values_prawn_imperfect <- function(specificity) {
  # Run the optimization using the log-likelihood function 'loglik_group_bb'
  optim_result <- optim(
    par = c(0, 0, 0),                 # Initial parameter guesses (log-scale)
    fn = loglik_group_bb,             # Log-likelihood function
    ty = summ.seafood.data$ty,        # Observed test results
    freq = summ.seafood.data$freq,    # Frequency of results
    b = 13,                           # Number of sampled groups
    m = 5,                            # Number of items per group
    M = 40,                           # Total items per consignment
    theta = Inf,                       # Overdispersion parameter to consider No clustering
    R = 1e4,                           # Number of static quantile from Beta distribution
    specificity = specificity,        # Input test specificity
    deviance = FALSE,                 # Do not return deviance
    hessian = FALSE,                  # No Hessian required
    control = list(
      factr = 1e-12,                  # High precision convergence
      fnscale = -1,                   # Maximize log-likelihood
      trace = TRUE                    # Print optimization progress
    ),
    method = "L-BFGS-B",              # Optimization method with bounds
    lower = c(log(0.0001), log(0.0001), log(0.10)),  # Parameter lower bounds
    upper = c(Inf, Inf, 0)            # Parameter upper bounds
  )

  # Extract and transform parameter estimates
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$sensitivity <- exp(optim_result$par[3])
  optim_result$mu <- optim_result$alpha / (optim_result$alpha + optim_result$beta)
  optim_result$specificity <- specificity

  # Return results as a data frame
  return(data.frame(
    alpha = optim_result$alpha,
    beta = optim_result$beta,
    mu.BB = optim_result$mu,
    sensitivity = optim_result$sensitivity,
    specificity = optim_result$specificity,
    value = optim_result$value,
    convergence = optim_result$convergence
  ))
}

# Define a range of specificity values to evaluate
specificity <- seq(0.9989, 1, by = 0.0000005)

# Apply the optimization safely with error handling
PLLF_Specificity_Free_Parameter_Prawn <- lapply(specificity, function(spcs) {
  tryCatch({
    store_specificity_values_prawn_imperfect(spcs)
  },
  error = function(e) {
    message("Error for specificity = ", spcs, ": ", e$message)
    return(NA)  # or NULL, depending on your downstream needs
  })
})

# Remove any failed results (NULLs) before combining
PLLF_Specificity_Free_Parameter_Prawn <-
  Filter(Negate(is.null), PLLF_Specificity_Free_Parameter_Prawn)

# Combine successful results into a single data frame
PLLF_Specificity_Free_Parameter_Prawn <-
  do.call(rbind, PLLF_Specificity_Free_Parameter_Prawn)

# Save the final combined results
save(
  PLLF_Specificity_Free_Parameter_Prawn,
  file = "PLLF_Specificity_Free_Parameter_Prawn_Re4.Rdata"
)



# ---------------------------------------------------------#
# PLLF for Sensitivity: Without restriction on Specificity
# ---------------------------------------------------------#

# MLE Estimates: alpha, beta and specificity
o.theta.inf.prawn.sp.sn <- optim(c(0,0,0),loglik_group_bb,ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,
                                 theta=Inf,m=5,M=40,sensitivity=1,deviance=FALSE,control=list(factr=1e-12,fnscale=-1,trace=T),
                                 R=1e4,hessian=FALSE,method = "L-BFGS-B",lower = c(log(0.0001),log(0.0001),log(0.9999)),upper = c(Inf,Inf,0))

# PLLF Calculation: Sensitivity  without restricting Specificity
# Function to run the optimization for a given Sensitivity value
store_sensitivity_values_prawn_imperfect <- function(sensitivity) {
  # Run the optimization using the log-likelihood function 'loglik_group_bb'
  optim_result <- optim(c(0, 0, 0), loglik_group_bb,
                        ty = summ.seafood.data$ty,
                        freq = summ.seafood.data$freq,
                        b = 13,
                        m = 5,
                        M = 40,
                        theta = Inf,
                        R = 1e4,
                        # R = 3e4,
                        # sensitivity = 0.8,
                        sensitivity = sensitivity,
                        deviance = FALSE,
                        hessian = FALSE,
                        control=list(factr=1e-12,fnscale=-1,trace=T),
                        method = "L-BFGS-B",
                        lower = c(log(0.0001),log(0.0001),log(0.99)),
                        upper = c(Inf,Inf,0)
  )

  # Store the parameter estimates
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$specificity <- exp(optim_result$par[3])

  optim_result$mu <- optim_result$alpha / (optim_result$alpha + optim_result$beta)
  optim_result$sensitivity <- sensitivity

  # Return the relevant values as a data frame
  return(data.frame(alpha = optim_result$alpha,
                    beta = optim_result$beta,
                    mu.BB = optim_result$mu,
                    sensitivity = optim_result$sensitivity,
                    specificity = optim_result$specificity,
                    #  mu.overall = optim_result$mu.overall,
                    value = optim_result$value,
                    convergence=optim_result$convergence
  ))
}

sensitivity  <- seq(0.45,1,by=0.0001)

PLLF_Sensitivity_Free_Parameter_Prawn <- lapply(sensitivity, function(sens) {
  tryCatch({
    store_sensitivity_values_prawn_imperfect(sens)
  }, error = function(e) {
    message("Error for sensitivity = ", sens, ": ", e$message)
    return(NA)  # or NULL, depending on your downstream needs
  })
})

# Remove any failed results (NULLs) before combining
PLLF_Sensitivity_Free_Parameter_Prawn <-
  Filter(Negate(is.null), PLLF_Sensitivity_Free_Parameter_Prawn)

# Combine successful results into a single data frame
PLLF_Sensitivity_Free_Parameter_Prawn <-
  do.call(rbind, PLLF_Sensitivity_Free_Parameter_Prawn)

# Save the final combined results
save(
  PLLF_Sensitivity_Free_Parameter_Prawn,
  file = "PLLF_Sensitivity_Free_Parameter_Prawn_Re4.Rdata"
)


#-------------------------------------------------------------
# PLLF: Sensitivity and specificity
# Supplementary Figure 2
#-------------------------------------------------------------

load("PLLF_Specificity_Free_Parameter_Prawn_Re4.Rdata")

PLLF_Specificity_Free_Parameter_Prawn <- PLLF_Specificity_Free_Parameter_Prawn[!is.na(PLLF_Specificity_Free_Parameter_Prawn$convergence) & PLLF_Specificity_Free_Parameter_Prawn$convergence==0,]
PLLF_Specificity_Free_Parameter_Prawn$logL <- PLLF_Specificity_Free_Parameter_Prawn$value

pllf.CI.estimate.specificity.95 <- pllf.CI.estimate(PLLF_Specificity_Free_Parameter_Prawn$value,parameter = PLLF_Specificity_Free_Parameter_Prawn$specificity,level = 0.95)
pllf.CI.estimate.specificity.90 <- pllf.CI.estimate(PLLF_Specificity_Free_Parameter_Prawn$value,parameter = PLLF_Specificity_Free_Parameter_Prawn$specificity,level = 0.90)

plot(PLLF_Specificity_Free_Parameter_Prawn$specificity,PLLF_Specificity_Free_Parameter_Prawn$value,typ="l")
plot(PLLF_Specificity_Free_Parameter_Prawn$sensitivity,PLLF_Specificity_Free_Parameter_Prawn$value,typ="l")

load("PLLF_Sensitivity_Free_Parameter_Prawn_Re4.Rdata")

PLLF_Sensitivity_Free_Parameter_Prawn <- PLLF_Sensitivity_Free_Parameter_Prawn[!is.na(PLLF_Sensitivity_Free_Parameter_Prawn$convergence) & PLLF_Sensitivity_Free_Parameter_Prawn$convergence==0,]

o.theta.inf.prawn.sn.sp.1 <- optim(c(0,0),loglik_group_bb,ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,
                                   theta=Inf,m=5,M=40,deviance=FALSE,control=list(reltol=1e-12,fnscale=-1),
                                   R=1e4,hessian=FALSE,sensitivity=1,specificity=1)

PLLF_Sensitivity_Free_Parameter_Prawn$value[PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity==1]<- o.theta.inf.prawn.sn.sp.1$value
PLLF_Sensitivity_Free_Parameter_Prawn$alpha[PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity==1]<- exp(o.theta.inf.prawn.sn.sp.1$par[1])
PLLF_Sensitivity_Free_Parameter_Prawn$beta[PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity==1]<- exp(o.theta.inf.prawn.sn.sp.1$par[2])
PLLF_Sensitivity_Free_Parameter_Prawn$mu.BB[PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity==1]<- exp(o.theta.inf.prawn.sn.sp.1$par[1])/(exp(o.theta.inf.prawn.sn.sp.1$par[1])+exp(o.theta.inf.prawn.sn.sp.1$par[2]))
PLLF_Sensitivity_Free_Parameter_Prawn$specificity[PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity==1]<- 1

PLLF_Sensitivity_Free_Parameter_Prawn<- PLLF_Sensitivity_Free_Parameter_Prawn[PLLF_Sensitivity_Free_Parameter_Prawn$specificity==1,]

PLLF_Sensitivity_Free_Parameter_Prawn$logL <- PLLF_Sensitivity_Free_Parameter_Prawn$value

pllf.CI.estimate.sensitivity.95 <- pllf.CI.estimate(PLLF_Sensitivity_Free_Parameter_Prawn$value,parameter = PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity,level = 0.95)
pllf.CI.estimate.sensitivity.90 <- pllf.CI.estimate(PLLF_Sensitivity_Free_Parameter_Prawn$value,parameter = PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity,level = 0.90)



# ---------------------------------------------#
# Supplementary Figure 2 : On Original scale logL
# ---------------------------------------------#

png(file="pllf CI sensitivity specifity free parameter.png", width = 10,height = 5,units = "in",res = 300)

par(mfrow=c(1,2), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 4, 4, 2))  # Adjusted margins

range_ll <- range(PLLF_Sensitivity_Free_Parameter_Prawn$logL)

vals.2 <- seq(round(range_ll[1],1), round(range_ll[2],2), length.out = 15)
formatted_labels_2 <- formatC(vals.2, format = "f", digits = 2)
vals.1 <- round(seq(min(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity),max(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity), length.out = 15),2)
formatted_labels_1 <- formatC(vals.1, format = "f", digits = 2)

plot(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity,PLLF_Sensitivity_Free_Parameter_Prawn$logL,typ="l",xlab="Sensitivity",ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1) # ,main="logL under perfect specificity"
axis(side = 2, at = vals.2, labels = formatted_labels_2 ,las = 1,cex.lab=0.6)
axis(side = 1, at = vals.1, labels = formatted_labels_1 , las = 2,cex.lab=0.6)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn$logL,na.rm=T),col=2,lty=3)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn$logL,na.rm=T)-3.84/2,col=4)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn$logL,na.rm=T)-2.71/2,col=5)
abline(v=pllf.CI.estimate.sensitivity.95$mle,col=2,lty=1)
abline(v=pllf.CI.estimate.sensitivity.95$mle.ll,col=4,lty=3)
abline(v=pllf.CI.estimate.sensitivity.95$mle.ul,col=4,lty=3)
abline(v=pllf.CI.estimate.sensitivity.90$mle.ll,col=5,lty=6)
abline(v=pllf.CI.estimate.sensitivity.90$mle.ul,col=5,lty=6)
# legend("bottomright",legend=c(paste("MLE=",mle.sensitivity),paste("LL=",ll.sensitivity),paste("UL=",ul.sensitivity)),bty = "n")


# Identify the x-value corresponding to the maximum y-value
max_x <- pllf.CI.estimate.sensitivity.95$mle
max_y <- max(PLLF_Sensitivity_Free_Parameter_Prawn$logL,na.rm=T)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
text(x=max_x-0.01, y=(max_y-1.5), labels=round(max_x, 5), pos=3, col=2, cex=0.7, srt = 90)
text(x=max_x+0.25, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.95$mle.ul, 5), pos=3, col=6, cex=0.5, srt = 90)
text(x=max_x-0.27, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.95$mle.ll, 5), pos=3, col=4, cex=0.5, srt = 90)
text(x=max_x-0.235, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.90$mle.ll, 5), pos=3, col=5, cex=0.5, srt = 90)

legend("bottomright",legend=c("0.95","0.90"),col=c(4,5),lty=c(3,6),bty="n",lwd=c(2,2))


range_ll <- range(PLLF_Specificity_Free_Parameter_Prawn$logL)
vals.2 <- round(seq(range_ll[1], range_ll[2], length.out = 15),2)
formatted_labels_2 <- formatC(vals.2, format = "f", digits = 2)
vals.1 <- round(seq(min(PLLF_Specificity_Free_Parameter_Prawn$specificity), max(PLLF_Specificity_Free_Parameter_Prawn$specificity), length.out = 15),5)
formatted_labels_1 <- formatC(vals.1, format = "f", digits = 4)

plot(PLLF_Specificity_Free_Parameter_Prawn$specificity,PLLF_Specificity_Free_Parameter_Prawn$logL,typ="l",xlab="Specificity",ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1) # ,main="logL under perfect specificity"
axis(side = 2, at = vals.2, labels = formatted_labels_2, las = 1,cex.lab=0.6)
axis(side = 1, at = vals.1, labels = formatted_labels_1, las = 2,cex.lab=0.6)
abline(h=max(PLLF_Specificity_Free_Parameter_Prawn$logL,na.rm=T),col=2,lty=3)
abline(h=max(PLLF_Specificity_Free_Parameter_Prawn$logL,na.rm=T)-3.84/2,col=4)
abline(h=max(PLLF_Specificity_Free_Parameter_Prawn$logL,na.rm=T)-2.71/2,col=5)
abline(v=pllf.CI.estimate.specificity.95$mle,col=2,lty=1)
abline(v=pllf.CI.estimate.specificity.95$mle.ll,col=4,lty=2)
abline(v=pllf.CI.estimate.specificity.95$mle.ul,col=4,lty=2)
abline(v=pllf.CI.estimate.specificity.90$mle.ll,col=5,lty=2)
# legend("bottomright",legend=c(paste("MLE=",mle.sensitivity),paste("LL=",ll.sensitivity),paste("UL=",ul.sensitivity)),bty = "n")


# Identify the x-value corresponding to the maximum y-value
max_x <- pllf.CI.estimate.specificity.95$mle
max_y <- max(PLLF_Specificity_Free_Parameter_Prawn$logL,na.rm=T)

# Display the x-value on the plot
text(x=max_x-0.00002, y=(max_y-15), labels=round(max_x, 5), pos=3, col=2, cex=0.7, srt = 90)
# text(x=max_x+0.25, y=(max_y-2), labels=round(pllf.CI.estimate.specificity.95$mle.ul, 5), pos=3, col=6, cex=0.9, srt = 90)
text(x=max_x-0.00015, y=(max_y-15), labels=round(pllf.CI.estimate.specificity.95$mle.ll, 5), pos=3, col=4, cex=0.5, srt = 90)
text(x=max_x-0.00011, y=(max_y-15), labels=round(pllf.CI.estimate.specificity.90$mle.ll, 5), pos=3, col=5, cex=0.5, srt = 90)

legend("bottomleft",legend=c("0.95","0.90"),col=c(4,5),lty=c(3,6),bty="n",lwd=c(2,2))

dev.off()

# ---------------------------------------------#
# Supplementary Figure 2 : On Re-scaled logL  (logL - maxlogL)
# ---------------------------------------------#



png(file="pllf CI sensitivity specifity free parameter Rescaled.png", width = 10,height = 5,units = "in",res = 300)

par(mfrow=c(1,2), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 4, 4, 2))  # Adjusted margins


PLLF_Sensitivity_Free_Parameter_Prawn$logL_rescaled <- PLLF_Sensitivity_Free_Parameter_Prawn$logL - max(PLLF_Sensitivity_Free_Parameter_Prawn$logL,na.rm=T)

range_ll <- range(PLLF_Sensitivity_Free_Parameter_Prawn$logL_rescaled)

step <- -0.25
ymax <- 0
ymin <- floor(range_ll[1]*10)/10  # round down to nearest 0.1
vals.2 <- seq(from = ymax, to = ymin, by = step)
formatted_labels_2 <- formatC(vals.2, format = "f", digits = 2)

step <- 0.05
xmax <- max(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity)
xmin <- min(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity) # round down to nearest 0.1
vals.1 <- seq(from = xmin, to = xmax, by = step)
formatted_labels_1 <- formatC(vals.1, format = "f", digits = 2)

plot(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity,
     PLLF_Sensitivity_Free_Parameter_Prawn$logL_rescaled,
     typ="l",
     xlab="Sensitivity",ylab="logL - max(logL)",
     yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1) # ,main="logL under perfect specificity"
axis(side = 2, at = vals.2, labels = formatted_labels_2 ,las = 1,cex.lab=0.6,tcl = -0.3,cex.axis=0.9)
axis(side = 1, at = vals.1, labels = formatted_labels_1 , las = 2,cex.lab=0.6,tcl = -0.3,cex.axis=0.9)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn$logL_rescaled,na.rm=T),col=2,lty=3)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn$logL_rescaled,na.rm=T)-3.84/2,col=4)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn$logL_rescaled,na.rm=T)-2.71/2,col=5)

legend("bottomright",legend=c("0.95","0.90"),col=c(4,5),lty=c(3,6),bty="n",lwd=c(2,2))

PLLF_Specificity_Free_Parameter_Prawn$logL_rescaled <- PLLF_Specificity_Free_Parameter_Prawn$logL -
                                                        max(PLLF_Specificity_Free_Parameter_Prawn$logL,na.rm=T)
range_ll <- range(PLLF_Specificity_Free_Parameter_Prawn$logL_rescaled)

step <- -2.5
ymax <- 0
ymin <- floor(range_ll[1]*10)/10  # round down to nearest 0.1
vals.2 <- seq(from = ymax, to = ymin, by = step)
formatted_labels_2 <- formatC(vals.2, format = "f", digits = 2)

step <- 8e-05
xmax <- max(PLLF_Specificity_Free_Parameter_Prawn$specificity)
xmin <- min(PLLF_Specificity_Free_Parameter_Prawn$specificity)
vals.1 <- seq(from = xmin, to = xmax, length.out = 15)
formatted_labels_1 <- formatC(vals.1, format = "f", digits = 4)

plot(PLLF_Specificity_Free_Parameter_Prawn$specificity,
     PLLF_Specificity_Free_Parameter_Prawn$logL_rescaled,
     typ="l",xlab="Specificity",
     ylab="logL - max(logL)",
     yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1) # ,main="logL under perfect specificity"
axis(side = 2, at = vals.2, labels = formatted_labels_2, las = 1,cex.lab=0.6,tcl = -0.3,cex.axis=0.9)
axis(side = 1, at = vals.1, labels = formatted_labels_1, las = 2,cex.lab=0.6,tcl = -0.3,cex.axis=0.9)
abline(h=max(PLLF_Specificity_Free_Parameter_Prawn$logL_rescaled,na.rm=T),col=2,lty=3)
abline(h=max(PLLF_Specificity_Free_Parameter_Prawn$logL_rescaled,na.rm=T)-3.84/2,col=4)
abline(h=max(PLLF_Specificity_Free_Parameter_Prawn$logL_rescaled,na.rm=T)-2.71/2,col=5)

legend("bottomright",legend=c("0.95","0.90"),col=c(4,5),lty=c(3,6),bty="n",lwd=c(2,2))

dev.off()

