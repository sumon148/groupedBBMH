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


# -------------------------------------------------------------@
# Figure 2: JASA PAper
# Profile likelihood for cutoff values
# Data: Deidentified Frozen Seafood Data
# Exact BB Model: MLE method and PLLF
# Using Generailzed log-likelihood function
# Target Parameters: PLLF of Minimum Positive Propensity (k)
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

# --------------------------------------------------------------
# PLLF of Minimum Positive Propensity: Perfect Sensitivity
# --------------------------------------------------------------


store_theta_values_prawn <- function(cut_off) {
  # Run the optimization
  optim_result <- optim(c(0, 0), loglik_group_trbb,
                        ty = summ.prawn.data$ty,
                        freq = summ.prawn.data$freq,
                        b = 13,
                        m = 5,
                        M = 40,
                        theta = Inf,
                        R = 1e5, # R = 3e4
                        sensitivity = 1,
                        specificity = 1,
                        cutoff = cut_off,
                        deviance = FALSE,
                        control = list(reltol = 1e-12, fnscale = -1),
                        hessian = FALSE)

  # Store the parameter estimates
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$cutoff <- cut_off
  optim_result$mu.TrBB <- optim_result$alpha / (optim_result$alpha + optim_result$beta)

  # Return the relevant values as a data frame
  return(data.frame(cutoff = optim_result$cutoff,
                    alpha = optim_result$alpha,
                    beta = optim_result$beta,
                    mu.BB = optim_result$mu.TrBB,
                    #  mu.overall = optim_result$mu.overall,
                    value = optim_result$value,
                    convergence=optim_result$convergence
  ))
}

cutoff_values <- seq(0.000, 0.03, by = 0.00001)

# Use lapply to apply the function to each pi value and store results as a list of dataframes
prawn_results_list <- lapply(cutoff_values, store_theta_values_prawn)

# Combine the list of dataframes into a single dataframe
prawn_results_df_TrBB <- do.call(rbind, prawn_results_list)

save(prawn_results_df_TrBB,file="prawn_results_df_TrBB_upto_3_percent_Re5.Rdata")


load("prawn_results_df_TrBB_upto_3_percent_Re5.Rdata")
pllf.CI.estimate.cutoff.99 <- pllf.CI.estimate(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB$value,level=0.99)
pllf.CI.estimate.cutoff.95 <- pllf.CI.estimate(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB$value,level=0.95)
pllf.CI.estimate.cutoff.90 <- pllf.CI.estimate(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB$value,level=0.90)
pllf.CI.estimate.cutoff.80 <- pllf.CI.estimate(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB$value,level=0.80)


png("pllf CI cutoff Prawn.png",width = 6,height = 6,units = "in",res=300)
par(mfrow=c(1,1), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 4, 4, 2))  # Adjusted margins

range_ll <- range(prawn_results_df_TrBB$value)
plot(prawn_results_df_TrBB$cutoff, prawn_results_df_TrBB$value,typ="l",
     ylab="logL",main="logL under perfect test",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1,
     xlab=expression(italic(k)))
axis(side = 2, at = seq(round(range_ll[1],1), round(range_ll[2],1), by = 1),las = 1)
axis(side = 1, at = seq(min(prawn_results_df_TrBB$cutoff), max(prawn_results_df_TrBB$cutoff), by = 0.0015),las = 3)
abline(h=max(prawn_results_df_TrBB$value,na.rm=T),col=2,lty=3)
# abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.99,1)/2,col=3,lty=3)
abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.95,1)/2,col=4,lty=3)
abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.90,1)/2,col=5,lty=6)
abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.80,1)/2,col=6,lty=6)
abline(v=pllf.CI.estimate.cutoff.95$mle,col=2,lty=1)
abline(v=pllf.CI.estimate.cutoff.95$mle.ul,col=3,lty=2)
abline(v=pllf.CI.estimate.cutoff.95$mle.ll,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.99$mle.ul,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.99$mle.ll,col=3,lty=2)
abline(v=pllf.CI.estimate.cutoff.90$mle.ul,col=4,lty=2)
abline(v=pllf.CI.estimate.cutoff.90$mle.ll,col=4,lty=2)
abline(v=pllf.CI.estimate.cutoff.80$mle.ul,col=6,lty=2)
abline(v=pllf.CI.estimate.cutoff.80$mle.ll,col=6,lty=2)

legend("bottomright",legend=c("0.95","0.90","0.80"),col=c(3,4,6),lty=c(2,2,2))


# Identify the x-value corresponding to the maximum y-value
max_x <- prawn_results_df_TrBB$cutoff[which.max(prawn_results_df_TrBB$value)]
max_y <- max(prawn_results_df_TrBB$value)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
text(x=max_x-0.0002, y=(max_y-40), labels=round(max_x, 4), pos=3, col=2, cex=0.9)
text(x=max_x+0.005, y=(max_y-20), labels=round(pllf.CI.estimate.cutoff.80$mle.ul, 4), pos=3, col=6, cex=0.9)
text(x=0.001, y=(max_y-10), labels=round(pllf.CI.estimate.cutoff.80$mle.ll, 4), pos=3, col=6, cex=0.9)
text(x=max_x+0.008, y=(max_y-30), labels=round(pllf.CI.estimate.cutoff.90$mle.ul, 4), pos=3, col=4, cex=0.9)
text(x=max_x+0.02, y=(max_y-40), labels=round(pllf.CI.estimate.cutoff.95$mle.ul, 4), pos=3, col=3, cex=0.9)

dev.off()

# --------------------------------------------------------------
# PLLF of Minimum Positive Propensity: Imperfect Sensitivity (0.8)
# --------------------------------------------------------------

store_theta_values_prawn_imperfect <- function(cut_off) {
  # Run the optimization
  optim_result <- optim(c(0, 0), loglik_group_trbb,
                        ty = summ.prawn.data$ty,
                        freq = summ.prawn.data$freq,
                        b = 13,
                        m = 5,
                        M = 40,
                        theta = Inf,
                        R = 1e5,
                        # R = 3e4,
                        sensitivity = 0.8,
                        specificity = 1,
                        cutoff = cut_off,
                        deviance = FALSE,
                        control = list(reltol = 1e-12, fnscale = -1),
                        hessian = FALSE)

  # Store the parameter estimates
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$cutoff <- cut_off
  optim_result$mu.TrBB <- optim_result$alpha / (optim_result$alpha + optim_result$beta)
  #  optim_result$mu.overall <- (1 - optim_result$pi) * optim_result$mu.BB

  # Return the relevant values as a data frame
  return(data.frame(cutoff = optim_result$cutoff,
                    alpha = optim_result$alpha,
                    beta = optim_result$beta,
                    mu.BB = optim_result$mu.TrBB,
                    #  mu.overall = optim_result$mu.overall,
                    value = optim_result$value,
                    convergence=optim_result$convergence
  ))
}

cutoff_values <- seq(0.000, 0.03, by = 0.00001)

# Use lapply to apply the function to each pi value and store results as a list of dataframes
prawn_results_list_imperfect <- lapply(cutoff_values, store_theta_values_prawn_imperfect)

# Combine the list of dataframes into a single dataframe
prawn_results_df_TrBB_imperfect <- do.call(rbind, prawn_results_list_imperfect)

save(prawn_results_df_TrBB_imperfect,file="prawn_results_df_TrBB_imperfect_upto_3_percent_Re5.Rdata")




# -----------------------------------------------------------
# Figure 2: On Original Scale of logL
# Under perfect sensitivity
# -----------------------------------------------------------

load("prawn_results_df_TrBB_upto_3_percent_Re5.Rdata")

pllf.CI.estimate.cutoff.99 <- pllf.CI.estimate(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB$value,level=0.99)
pllf.CI.estimate.cutoff.95 <- pllf.CI.estimate(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB$value,level=0.95)
pllf.CI.estimate.cutoff.90 <- pllf.CI.estimate(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB$value,level=0.90)
pllf.CI.estimate.cutoff.80 <- pllf.CI.estimate(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB$value,level=0.80)


png("pllf CI cutoff Prawn.png",width = 6,height = 6,units = "in",res=300)
 par(mfrow=c(1,1), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 4, 4, 2))  # Adjusted margins

#par(mfrow = c(1,1), mgp   = c(3, 0.5, 0),      # axis title and label spacing
#  oma   = c(0,0,0,0),
#  mar   = c(4, 4, 1, 1)      # reduced top (1) and right (1) margins
#)

range_ll <- range(prawn_results_df_TrBB$value)
plot(prawn_results_df_TrBB$cutoff, prawn_results_df_TrBB$value,typ="l",
     ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1,
     xlab=expression(italic(k))) # main="logL under perfect test",
axis(side = 2, at = seq(round(range_ll[1],1), round(range_ll[2],1), by = 1),las = 1)
axis(side = 1, at = seq(min(prawn_results_df_TrBB$cutoff), max(prawn_results_df_TrBB$cutoff), by = 0.0015),las = 3)
abline(h=max(prawn_results_df_TrBB$value,na.rm=T),col=2,lty=3)
# abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.99,1)/2,col=3,lty=3)
abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.95,1)/2,col=4,lty=3)
abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.90,1)/2,col=5,lty=6)
abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.80,1)/2,col=6,lty=6)
abline(v=pllf.CI.estimate.cutoff.95$mle,col=2,lty=1)
abline(v=pllf.CI.estimate.cutoff.95$mle.ul,col=3,lty=2)
abline(v=pllf.CI.estimate.cutoff.95$mle.ll,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.99$mle.ul,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.99$mle.ll,col=3,lty=2)
abline(v=pllf.CI.estimate.cutoff.90$mle.ul,col=4,lty=2)
abline(v=pllf.CI.estimate.cutoff.90$mle.ll,col=4,lty=2)
abline(v=pllf.CI.estimate.cutoff.80$mle.ul,col=6,lty=2)
abline(v=pllf.CI.estimate.cutoff.80$mle.ll,col=6,lty=2)

legend("bottomright",legend=c("0.95","0.90","0.80"),col=c(3,4,6),lty=c(2,2,2),bty = "n")


# Identify the x-value corresponding to the maximum y-value
max_x <- prawn_results_df_TrBB$cutoff[which.max(prawn_results_df_TrBB$value)]
max_y <- max(prawn_results_df_TrBB$value)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
text(x=max_x-0.0002, y=(max_y-2), labels=round(max_x, 5), pos=3, col=2, cex=0.9, srt = 90)
text(x=max_x+0.008, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.80$mle.ul, 5), pos=3, col=6, cex=0.9, srt = 90)
text(x=0.002, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.80$mle.ll, 5), pos=3, col=6, cex=0.9, srt = 90)
text(x=max_x-0.009, y=(max_y-3), labels=round(pllf.CI.estimate.cutoff.90$mle.ll, 5), pos=3, col=4, cex=0.9, srt = 90)
text(x=max_x+0.011, y=(max_y-3), labels=round(pllf.CI.estimate.cutoff.90$mle.ul, 5), pos=3, col=4, cex=0.9, srt = 90)
text(x=max_x+0.015, y=(max_y-3), labels=round(pllf.CI.estimate.cutoff.95$mle.ul, 5), pos=3, col=3, cex=0.9, srt = 90)

dev.off()


# -----------------------------------------------------------
# Figure 2: On Original Scale of logL
# Under imperfect sensitivity
# -----------------------------------------------------------

load("prawn_results_df_TrBB_imperfect_upto_3_percent_Re5.Rdata")
pllf.CI.estimate.cutoff.99.imperfect <- pllf.CI.estimate(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$value,level=0.99)
pllf.CI.estimate.cutoff.95.imperfect <- pllf.CI.estimate(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$value,level=0.95)
pllf.CI.estimate.cutoff.90.imperfect <- pllf.CI.estimate(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$value,level=0.90)
pllf.CI.estimate.cutoff.80.imperfect <- pllf.CI.estimate(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$value,level=0.80)

dim(prawn_results_df_TrBB_imperfect)

# LR test ------------------------#


logL_k0 <- prawn_results_df_TrBB_imperfect$value[prawn_results_df_TrBB_imperfect$cutoff==0]
LRT <- 2 * (prawn_results_df_TrBB_imperfect$value - logL_k0)
pvals <- 1 - pchisq(LRT, df = 1)  # No boundary correction
pvals_exact <- 0.5 * (1 - pchisq(LRT, df = 1)) # Boundary correction

png("LRT k equal zero.png",width = 5,height = 5,units = "in",res=300)
par(mfrow=c(1,1), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 4, 4, 2))  # Adjusted margins

plot(prawn_results_df_TrBB_imperfect$cutoff, pvals,
     type = "l", col = "blue", lwd = 2,
     ylim = c(0, 1),
     xlab = "k (cutoff)", ylab = "p-value",
     main = expression("LRT p-values for H"[0]*": k = 0"))

lines(prawn_results_df_TrBB_imperfect$cutoff, pvals_exact,
      col = "red", lwd = 2, lty = 2)

abline(h = 0.05, col = "darkgreen", lty = 3,lwd=2)  # 5% threshold
abline(h = 0.10, col = "pink", lty = 4,lwd=2)  # 10% threshold

legend("top",
       legend = c(
         expression(chi^2*"(1) p-value"),
         expression("corrected p-value")
       ),
       col = c("blue", "red"), lty = c(1,2), lwd = 2, bty="n",cex=0.6)

dev.off()


png("pllf CI cutoff Prawn imperfect.png",width = 5,height = 5,units = "in",res=300)
par(mfrow=c(1,1), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 4, 4, 2))  # Adjusted margins

range_ll <- range(prawn_results_df_TrBB_imperfect$value)

plot(prawn_results_df_TrBB_imperfect$cutoff, prawn_results_df_TrBB_imperfect$value,typ="l",
     ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1,
     xlab=expression(k)) # main="logL under imperfect test",
axis(side = 2, at = seq(round(range_ll[1],1), round(range_ll[2],1), by = 1),las = 1)
axis(side = 1, at = seq(min(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB_imperfect$cutoff), max(prawn_results_df_TrBB$cutoff,prawn_results_df_TrBB_imperfect$cutoff), by = 0.0015),las = 3)
abline(h=max(prawn_results_df_TrBB_imperfect$value,na.rm=T),col=2,lty=3)
# abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.99,1)/2,col=3,lty=3)
abline(h=max(prawn_results_df_TrBB_imperfect$value,na.rm=T)-qchisq(0.95,1)/2,col=4,lty=3)
abline(h=max(prawn_results_df_TrBB_imperfect$value,na.rm=T)-qchisq(0.90,1)/2,col=5,lty=6)
abline(h=max(prawn_results_df_TrBB_imperfect$value,na.rm=T)-qchisq(0.80,1)/2,col=6,lty=6)
abline(v=pllf.CI.estimate.cutoff.95.imperfect$mle,col=2,lty=1)
abline(v=pllf.CI.estimate.cutoff.95.imperfect$mle.ul,col=3,lty=2)
abline(v=pllf.CI.estimate.cutoff.95.imperfect$mle.ll,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.99$mle.ul,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.99$mle.ll,col=3,lty=2)
abline(v=pllf.CI.estimate.cutoff.90.imperfect$mle.ul,col=4,lty=2)
abline(v=pllf.CI.estimate.cutoff.90.imperfect$mle.ll,col=4,lty=2)
abline(v=pllf.CI.estimate.cutoff.80.imperfect$mle.ul,col=6,lty=2)
abline(v=pllf.CI.estimate.cutoff.80.imperfect$mle.ll,col=6,lty=2)

legend("topright",legend=c("0.95","0.90","0.80"),col=c(3,4,6),lty=c(2,2,2),bty="n",lwd=c(2,2,2))

# Identify the x-value corresponding to the maximum y-value
max_x <- prawn_results_df_TrBB_imperfect$cutoff[which.max(prawn_results_df_TrBB_imperfect$value)]
max_y <- max(prawn_results_df_TrBB_imperfect$value)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
text(x=max_x-0.0002, y=(max_y-2), labels=round(max_x, 5), pos=3, col=2, cex=0.9, srt = 90)
text(x=max_x+0.0105, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.80.imperfect$mle.ul, 5), pos=3, col=6, cex=0.9, srt = 90)
text(x=max_x-0.008, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.80.imperfect$mle.ll, 5), pos=3, col=6, cex=0.9, srt = 90)
text(x=max_x+0.014, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.90.imperfect$mle.ul, 5), pos=3, col=4, cex=0.9, srt = 90)
text(x=max_x+0.017, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.95.imperfect$mle.ul, 5), pos=3, col=3, cex=0.9, srt = 90)


dev.off()


# -----------------------------------------------------------
# Figure 2: On Re-Scaled logL: logL - max(logL)
# Under imperfect sensitivity
# -----------------------------------------------------------

load("prawn_results_df_TrBB_imperfect_upto_3_percent_Re5.Rdata")
pllf.CI.estimate.cutoff.99.imperfect <- pllf.CI.estimate(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$value,level=0.99)
pllf.CI.estimate.cutoff.95.imperfect <- pllf.CI.estimate(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$value,level=0.95)
pllf.CI.estimate.cutoff.90.imperfect <- pllf.CI.estimate(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$value,level=0.90)
pllf.CI.estimate.cutoff.80.imperfect <- pllf.CI.estimate(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$value,level=0.80)


prawn_results_df_TrBB_imperfect$value_rescaled <- prawn_results_df_TrBB_imperfect$value - max(prawn_results_df_TrBB_imperfect$value,na.rm=TRUE)

png("pllf CI cutoff Prawn imperfect rescaled.png",width = 6,height = 5,units = "in",res=300)
# par(mfrow=c(1,1), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 5.5, 2, 2))  # Adjusted margins
par(mfrow = c(1, 1), mgp = c(4, 1, 0), oma = c(0, 0, 0, 0), mar = c(5, 5.5, 2, 2))
range_ll <- range(prawn_results_df_TrBB_imperfect$value_rescaled)

step <- -0.25
ymax <- 0
ymin <- floor(range_ll[1]*10)/10  # round down to nearest 0.1
y_ticks <- seq(from = ymax, to = ymin, by = step)

plot(prawn_results_df_TrBB_imperfect$cutoff, prawn_results_df_TrBB_imperfect$value_rescaled,typ="l",
     ylab="logL - max(logL)",yaxt="n",xaxt="n",cex.lab=1.5,cex.axis=0.9,
     xlab=expression(italic(k)) ) # main="logL under imperfect test",

axis(side = 2, at = y_ticks,las = 1, cex.axis = 0.9,tcl = -0.5)
axis(side = 1, at = seq(min(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$cutoff), max(prawn_results_df_TrBB_imperfect$cutoff,prawn_results_df_TrBB_imperfect$cutoff), by = 0.0015),
     las = 3,tcl = -0.5,
     cex.axis = 0.9)
abline(h=max(prawn_results_df_TrBB_imperfect$value_rescaled,na.rm=T),col=2,lty=2)

# abline(h=max(prawn_results_df_TrBB$value,na.rm=T)-qchisq(0.99,1)/2,col=3,lty=3)
# abline(h=max(prawn_results_df_TrBB_imperfect$value_rescaled,na.rm=T)-qchisq(0.95,1)/2,col=4,lty=3)
# abline(h=max(prawn_results_df_TrBB_imperfect$value_rescaled,na.rm=T)-qchisq(0.90,1)/2,col=5,lty=6)
# abline(h=max(prawn_results_df_TrBB_imperfect$value_rescaled,na.rm=T)-qchisq(0.80,1)/2,col=6,lty=6)
# abline(v=pllf.CI.estimate.cutoff.95.imperfect$mle,col=2,lty=1)
#abline(v=pllf.CI.estimate.cutoff.95.imperfect$mle.ul,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.95.imperfect$mle.ll,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.99$mle.ul,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.99$mle.ll,col=3,lty=2)
#abline(v=pllf.CI.estimate.cutoff.90.imperfect$mle.ul,col=4,lty=2)
#abline(v=pllf.CI.estimate.cutoff.90.imperfect$mle.ll,col=4,lty=2)
# abline(v=pllf.CI.estimate.cutoff.80.imperfect$mle.ul,col=6,lty=2)
# abline(v=pllf.CI.estimate.cutoff.80.imperfect$mle.ll,col=6,lty=2)

# legend("topright",legend=c("0.95","0.90","0.80"),col=c(3,4,6),lty=c(2,2,2),bty="n",lwd=c(2,2,2))
# legend("bottomleft",legend=c("0.95","0.90"),col=c(4,5),lty=c(3,6),bty="n",lwd=c(2,2))

# Identify the x-value corresponding to the maximum y-value
# max_x <- prawn_results_df_TrBB_imperfect$cutoff[which.max(prawn_results_df_TrBB_imperfect$value_rescaled)]
# max_y <- max(prawn_results_df_TrBB_imperfect$value_rescaled)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
# text(x=max_x-0.0002, y=(max_y-2), labels=round(max_x, 5), pos=3, col=2, cex=0.9, srt = 90)
# text(x=max_x+0.0105, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.80.imperfect$mle.ul, 5), pos=3, col=6, cex=0.9, srt = 90)
# text(x=max_x-0.008, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.80.imperfect$mle.ll, 5), pos=3, col=6, cex=0.9, srt = 90)
# text(x=max_x+0.014, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.90.imperfect$mle.ul, 5), pos=3, col=4, cex=0.9, srt = 90)
# text(x=max_x+0.017, y=(max_y-2), labels=round(pllf.CI.estimate.cutoff.95.imperfect$mle.ul, 5), pos=3, col=3, cex=0.9, srt = 90)


dev.off()



