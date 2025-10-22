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
# Compiling Results for Figure 3 & 4: JASA PAper
# Posterior Predictive Check and Dunn-Smyth Diagnostic
# Using Prawn data
# Exact BB Model with \Delta=0.80: MH in Paper
# Exact BB Model with \Delta=1.0: MH in Supplementary File
# Approximate BB Model with \Delta=1.0: HMC in Supplementary File
# Relevant Supplementary Table will be from here
#-------------------------------------------------------------@

library(ggplot2)
library(ggpubr)


ppc_barplot_fixed_y_scale <- function(ty,yrep,trial){
  # ToDo: Convert into one plot and put legend in the plot
  # ty= Vector of counts
  # yrep = Draws of counts from posterior distribution.
  # trial = Number of trials


  iter = dim(yrep)[1]
  tn = length(ty) # tn = Realization size
  yrep.mutate<-matrix(NA,iter,(trial+1))
  colnames(yrep.mutate) <- as.character(c(0:trial))
  list.yrep <- apply(yrep,1,function(x){data.frame(table(x))})


  for (i in 1:iter){
    vector1 <- c(0:trial)
    position <- vector1%in%list.yrep[[i]]$x
    yrep.mutate[i,position] <- list.yrep[[i]]$Freq
  }

  yrep.mutate[is.na(yrep.mutate)] <- 0

  y.counts <- data.frame((table(ty)))

  mean.y <- apply(yrep.mutate,2,function(x){mean(x,na.rm=T)})
  median.y <- apply(yrep.mutate,2,function(x){median(x,na.rm=T)})
  sd.y <- apply(yrep.mutate,2,function(x){sd(x,na.rm=T)})
  y.025 <- apply(yrep.mutate,2,function(x){quantile(x,0.025,na.rm=T)})
  y.975 <- apply(yrep.mutate,2,function(x){quantile(x,0.975,na.rm=T)})


  count.summary <- data.frame(x=as.character(c(0:trial)),mean.y=mean.y,median.y=median.y,sd.y=sd.y,y.025=y.025,y.975=y.975)

  m <- (count.summary$x%in%y.counts$ty)
  count.summary$y.observed <- NA
  count.summary$y.observed[m] <-  y.counts$Freq

  dodge <- position_dodge(width=1)
  p1 <- ggplot(count.summary[count.summary$x=="0",], aes(x=x, y=y.observed))
  p1 <- p1 +
    geom_col(position = "dodge",fill="skyblue",just = 0.5,width = 0.25) +
    geom_errorbar(aes(ymin = y.025, ymax = y.975,fill=NULL), position = "dodge", width = 0.05)+
    geom_point(aes(x=x,y=mean.y),size=1.5)+
    #    geom_point(aes(x=x,y=median.y),size=2,shape=8)+
    scale_x_discrete("Infected bags", limits = as.character(c(0)))+ylab("Frequency") #+ggtitle("Counts with 95% CI")

  p2 <- ggplot(count.summary[count.summary$x%in%as.character(1:trial),], aes(x=x, y=y.observed))
  p2 <- p2 +
    geom_col(position = "dodge",fill="skyblue") +
    geom_errorbar(aes(ymin = y.025, ymax = y.975,fill=NULL), position = "dodge", width = 0.25)+
    geom_point(aes(x=x,y=mean.y),size=1.5)+
    #    geom_point(aes(x=x,y=median.y),size=2,shape=8)+
    scale_x_discrete("Infected bags", limits = as.character(c(1:trial)))+ylab("Frequency") # +ggtitle("Counts with 95% CI")

  p.counts <-  ggarrange(p1,p2,nrow=1)
  p.counts <- annotate_figure(p.counts, top = text_grob(paste("Counts (out of ",tn,") with 95% CI"), color = "red", face = "bold", size = 14))


  # Extra -----------#
  count.summary.dup <- count.summary
  count.summary.dup$mean.y1 <- NA
  count.summary.dup$mean.y1[count.summary.dup$x=="0"] <- count.summary.dup$mean.y[count.summary.dup$x=="0"]
  count.summary.dup$mean.y2 <- NA
  count.summary.dup$mean.y2[count.summary.dup$x!="0"] <- count.summary.dup$mean.y[count.summary.dup$x!="0"]
  count.summary.dup$y1.observed <- NA
  count.summary.dup$y2.observed <- NA
  count.summary.dup$y1.observed[count.summary.dup$x=="0"]  <- count.summary.dup$y.observed[count.summary.dup$x=="0"]
  count.summary.dup$y2.observed[count.summary.dup$x!="0"] <- count.summary.dup$y.observed[count.summary.dup$x!="0"]
  count.summary.dup$y1.025 <- NA
  count.summary.dup$y2.025 <- NA
  count.summary.dup$y1.025[count.summary.dup$x=="0"]  <- count.summary.dup$y.025[count.summary.dup$x=="0"]
  count.summary.dup$y2.025[count.summary.dup$x!="0"] <- count.summary.dup$y.025[count.summary.dup$x!="0"]
  count.summary.dup$y1.975 <- NA
  count.summary.dup$y2.975 <- NA
  count.summary.dup$y1.975[count.summary.dup$x=="0"]  <- count.summary.dup$y.975[count.summary.dup$x=="0"]
  count.summary.dup$y2.975[count.summary.dup$x!="0"] <- count.summary.dup$y.975[count.summary.dup$x!="0"]

  # Create the plot with geom_col for both variables
  count.summary.dup$x <- factor(count.summary.dup$x,levels = as.character(c(0:trial)))

  ratio.zero.nonzero <- sum(count.summary$y.observed[count.summary$x!="0"],na.rm=T)/count.summary$y.observed[count.summary$x=="0"]

  # ggplot(count.summary.dup, aes(x = x)) +
  # geom_col(aes(y = y1.observed), fill = "blue", alpha = 0.5, width = 0.4) +  # First y variable (y1) as bars
  # geom_col(aes(y = y2.observed*1/0.007), fill = "red", alpha = 0.5, width = 0.4) +  # Second y variable (y2), scaled and bars
  # scale_y_continuous(
  #   name = "Number of zero counts",         # Primary y-axis
  #   sec.axis = sec_axis(~.*0.007, name = "Number of non-zero counts")  # Secondary y-axis
  # ) +
  # theme_minimal() +
  # theme(
  #   axis.title.y.left = element_text(color = "blue"),
  #   axis.title.y.right = element_text(color = "red"),
  #   legend.position = "none"
  # ) +
  # labs(x = "Number of positive test", title = "Number of zero and non-zero counts")


  p.zero.non.zero <- ggplot(count.summary.dup, aes(x = x)) +
    geom_col(aes(y = y1.observed), fill = "blue", alpha = 0.5, width = 0.4) +  # First y variable (y1) as bars
    geom_col(aes(y = y2.observed*3000/25), fill = "red", alpha = 0.5, width = 0.4) +  # Second y variable (y2), scaled and bars
    # Error bars for y1.observed (zero counts)
    geom_errorbar(aes(ymin = y1.025, ymax = y1.975), width = 0.2, color = "blue") +
    # Error bars for y2.observed (non-zero counts) - scaled
    geom_errorbar(aes(ymin = y2.025 * 3000/25, ymax = y2.975 * 3000/25), width = 0.2, color = "red") +

    # For zero counts posterior mean
    geom_point(data = count.summary.dup[x==0,], aes(x = x, y = mean.y),
               color = "blue", size = 1, shape = 18) +

    # For non-zero counts posterior mean (scaled)
    geom_point(data = count.summary.dup[x>0,], aes(x = x, y = mean.y * 3000 / 25),
               color = "red", size = 1, shape = 18)  +

    scale_y_continuous(
      name = "Number of zero counts",             # Primary y-axis
      limits = c(0, 3000),                        # Fixed limit for primary y-axis
      breaks = seq(0, 3000, 500),                 # Breaks for primary y-axis
      sec.axis = sec_axis(
        ~ . * 25 / 3000,                          # Transformation for secondary axis
        name = "Number of non-zero counts",
        breaks = seq(0, 25, 5)                    # Breaks for secondary y-axis
      )
    ) +
    theme_minimal() +
    theme(
      axis.title.y.left = element_text(color = "blue"),
      axis.title.y.right = element_text(color = "red"),
      legend.position = "none"
    ) +
    labs(x = "Infected groups", title = paste("Counts (out of ",tn,") with 95% CI"))


  # ------------------------#


  mean.prop <- apply(yrep.mutate/tn,2,function(x){mean(x,na.rm=T)})
  median.prop <- apply(yrep.mutate/tn,2,function(x){median(x,na.rm=T)})
  sd.prop <- apply(yrep.mutate/tn,2,function(x){sd(x,na.rm=T)})
  prop.025 <- apply(yrep.mutate/tn,2,function(x){quantile(x,0.025,na.rm=T)})
  prop.975 <- apply(yrep.mutate/tn,2,function(x){quantile(x,0.975,na.rm=T)})

  y.prop <- data.frame(prop.table(table(ty)))
  colnames(y.prop) <- c("x","Prop")
  prop.summary <- data.frame(x=as.character(c(0:trial)),mean.prop=mean.prop,median.prop=median.prop,sd.prop=sd.prop,prop.025=prop.025,prop.975=prop.975)

  m <- (prop.summary$x%in%y.prop$x)
  prop.summary$y.observed <- NA
  prop.summary$y.observed[m] <-  y.prop$Prop

  dodge <- position_dodge(width=1)
  p1 <- ggplot(prop.summary[prop.summary$x=="0",], aes(x=x, y=y.observed))
  p1 <- p1 +
    geom_col(position = "dodge",fill="skyblue",just = 0.5,width = 0.25) +
    geom_errorbar(aes(ymin = prop.025, ymax = prop.975,fill=NULL), position = "dodge", width = 0.05)+
    geom_point(aes(x=x,y=mean.prop),size=1.5)+
    #    geom_point(aes(x=x,y=median.y),size=2,shape=8)+
    scale_x_discrete("Infected bags", limits = as.character(c(0)))+ylab("Proportion")+ylim(0,1.05) # +ggtitle("Proportion with 95% CI")

  p2 <- ggplot(prop.summary[prop.summary$x%in%as.character(1:trial),], aes(x=x, y=y.observed))
  p2 <- p2 +
    geom_col(position = "dodge",fill="skyblue") +
    geom_errorbar(aes(ymin = prop.025, ymax = prop.975,fill=NULL), position = "dodge", width = 0.25)+
    geom_point(aes(x=x,y=mean.prop),size=1.5)+
    #    geom_point(aes(x=x,y=median.y),size=2,shape=8)+
    scale_x_discrete("Infected bags", limits = as.character(c(1:trial)))+ylab("Proportion") # +ggtitle("Proportion with 95% CI")


  p.prop <-  ggarrange(p1,p2,nrow=1)
  p.prop <- annotate_figure(p.prop, top = text_grob("Proportion with 95% CI", color = "red", face = "bold", size = 14))


  # Extra -----------#
  prop.summary.dup <- prop.summary
  prop.summary.dup$mean.y1 <- NA
  prop.summary.dup$mean.y1[prop.summary.dup$x=="0"] <- prop.summary.dup$mean.prop[prop.summary.dup$x=="0"]
  prop.summary.dup$mean.y2 <- NA
  prop.summary.dup$mean.y2[prop.summary.dup$x!="0"] <- prop.summary.dup$mean.prop[prop.summary.dup$x!="0"]
  prop.summary.dup$y1.observed <- NA
  prop.summary.dup$y2.observed <- NA
  prop.summary.dup$y1.observed[prop.summary.dup$x=="0"]  <- prop.summary.dup$y.observed[prop.summary.dup$x=="0"]
  prop.summary.dup$y2.observed[prop.summary.dup$x!="0"] <- prop.summary.dup$y.observed[prop.summary.dup$x!="0"]
  prop.summary.dup$y1.025 <- NA
  prop.summary.dup$y2.025 <- NA
  prop.summary.dup$y1.025[prop.summary.dup$x=="0"]  <- prop.summary.dup$prop.025[prop.summary.dup$x=="0"]
  prop.summary.dup$y2.025[prop.summary.dup$x!="0"] <- prop.summary.dup$prop.025[prop.summary.dup$x!="0"]
  prop.summary.dup$y1.975 <- NA
  prop.summary.dup$y2.975 <- NA
  prop.summary.dup$y1.975[prop.summary.dup$x=="0"]  <- prop.summary.dup$prop.975[prop.summary.dup$x=="0"]
  prop.summary.dup$y2.975[prop.summary.dup$x!="0"] <- prop.summary.dup$prop.975[prop.summary.dup$x!="0"]

  # Create the plot with geom_col for both variables
  prop.summary.dup$x <- factor(prop.summary.dup$x,levels = as.character(c(0:trial)))

  ratio.zero.nonzero.prop <- (1-prop.summary$y.observed[count.summary$x=="0"])/(prop.summary$y.observed[count.summary$x=="0"])


  prop.zero.non.zero <- ggplot(prop.summary.dup, aes(x = x)) +
    geom_col(aes(y = y1.observed), fill = "blue", alpha = 0.5, width = 0.4) +  # First y variable (y1) as bars
    geom_col(aes(y = y2.observed*1/ratio.zero.nonzero.prop), fill = "red", alpha = 0.5, width = 0.4) +  # Second y variable (y2), scaled and bars

    # Error bars for y1.observed (zero counts)
    geom_errorbar(aes(ymin = y1.025, ymax = y1.975), width = 0.2, color = "blue") +
    # Error bars for y2.observed (non-zero counts) - scaled
    geom_errorbar(aes(ymin = y2.025 * 1/ratio.zero.nonzero.prop, ymax = y2.975 * 1/ratio.zero.nonzero.prop), width = 0.2, color = "red") +

    scale_y_continuous(
      name = "Proportion of zero counts",         # Primary y-axis
      sec.axis = sec_axis(~.*ratio.zero.nonzero.prop, name = "Proportion of non-zero counts")  # Secondary y-axis
    ) +
    theme_minimal() +
    theme(
      axis.title.y.left = element_text(color = "blue"),
      axis.title.y.right = element_text(color = "red"),
      legend.position = "none"
    ) +
    labs(x = "Infected groups", title = paste("Proportion with 95% CI"))



  # Zero and n-zero counts

  colnames(yrep.mutate)
  yrep.mutate.zero <- yrep.mutate[,"0"]
  yrep.mutate.non.zero <- apply(yrep.mutate[,colnames(yrep.mutate)!="0"],1,sum)
  yrep.mutate.zero.non.zero <- data.frame(zero=yrep.mutate.zero,non.zero=yrep.mutate.non.zero)
  mean.y.zero <- apply(yrep.mutate.zero.non.zero,2,function(x){mean(x,na.rm=T)})
  median.y.zero <- apply(yrep.mutate.zero.non.zero,2,function(x){median(x,na.rm=T)})
  sd.y.zero <- apply(yrep.mutate.zero.non.zero,2,function(x){sd(x,na.rm=T)})
  y.025.zero <- apply(yrep.mutate.zero.non.zero,2,function(x){quantile(x,0.025,na.rm=T)})
  y.975.zero <- apply(yrep.mutate.zero.non.zero,2,function(x){quantile(x,0.975,na.rm=T)})
  y.observed.zero <- c(y.counts$Freq[y.counts$ty=="0"],sum(y.counts$Freq[y.counts$ty!="0"]))
  count.zero.summary <- data.frame(x=c("Zero","Non-zero"),mean.y=mean.y.zero,median.y=median.y.zero,sd.y=sd.y.zero,y.025=y.025.zero,y.975=y.975.zero)

  count.zero.summary$y.observed <- y.observed.zero


  p1 <- ggplot(count.zero.summary, aes(x=x, y=y.observed))
  p1 <- p1 +
    geom_col(position = "dodge",fill="skyblue") +
    geom_errorbar(aes(ymin = y.025, ymax = y.975,fill=NULL), position = "dodge", width = 0.25)+
    geom_point(aes(x=x,y=mean.y),size=1.5)+
    #    geom_point(aes(x=x,y=median.y),size=2,shape=8)+
    scale_x_discrete("Infected bags",limits = c("Zero","Non-zero"))+ylab("Frequency") +ggtitle("Zero and non-zero counts with 95% CI")

  count.zero.proportion <- count.zero.summary
  count.zero.proportion$mean.y <- count.zero.summary$mean.y/sum(count.zero.summary$y.observed)
  count.zero.proportion$median.y.zero <- count.zero.summary$median.y/sum(count.zero.summary$y.observed)
  count.zero.proportion$y.025.zero <- count.zero.summary$y.025/sum(count.zero.summary$y.observed)
  count.zero.proportion$y.975.zero <- count.zero.summary$y.975/sum(count.zero.summary$y.observed)
  count.zero.proportion$y.prop <- count.zero.summary$y.observed/sum(count.zero.summary$y.observed)

  p2 <- ggplot(count.zero.proportion, aes(x=x, y=y.prop))
  p2 <- p2 +
    geom_col(position = "dodge",fill="skyblue") +
    geom_errorbar(aes(ymin = y.025.zero, ymax = y.975.zero,fill=NULL), position = "dodge", width = 0.25)+
    geom_point(aes(x=x,y=mean.y),size=1.5)+
    #    geom_point(aes(x=x,y=median.y),size=2,shape=8)+
    scale_x_discrete("Infected bags",limits = c("Zero","Non-zero"))+ylab("Proportion") +ggtitle("Zero and non-zero proportions with 95% CI")


  p.prop.zero <-  ggarrange(p1,p2,nrow=1)

  list(count.summary=count.summary,prop.summary=prop.summary,p.counts=p.counts,p.prop=p.prop,p.prop.zero=p.prop.zero,count.zero.summary=count.zero.summary,p.zero.non.zero=p.zero.non.zero,prop.zero.non.zero=prop.zero.non.zero)
}

#-------------------------------------------------------------
# Prawn Data
#-------------------------------------------------------------@

df.prawn <- read.csv("JASA Submission/deidentified_frozen_seafood.csv")
df.prawn$Yfac <- factor(df.prawn$numPos,levels=c(0:13))
x <- as.numeric(names(table(df.prawn$Yfac)))
freq <- as.numeric(table(df.prawn$Yfac))
size <- 13 # Number of bags (b)
prawn.data <- data.frame(ty=df.prawn$Yfac,n=rep(size,length(df.prawn$Yfac)))
prawn.data$ID <- c(1:dim(prawn.data)[1])
prawn.data$ty <- as.numeric(paste(prawn.data$ty))
prawn.data <- prawn.data[order(prawn.data$ty),]

summ.prawn.data <- as.data.frame(table(prawn.data$ty))
colnames(summ.prawn.data) <- c("ty","freq")
summ.prawn.data$ty <- as.numeric(paste(summ.prawn.data$ty))

# ----------------------------------------------------------------------------------
# Posterior Predictive Check: Imperfect Sensitivity (Sensitivity = 0.80)
# Model: Exact Model using Metropolis-Hastings Algorithm with known Delta = 0.8
# Data: Prawn Infection Dataset
# Purpose: Generate posterior predictive samples under imperfect test sensitivity
# ----------------------------------------------------------------------------------

# Load posterior samples from MCMC output (Metropolis-Hastings)
load("JASA Submission/MH.alpha.mu.sigma.20.known.sn.80.Prawn.Rdata")

# Number of test groups (D) from observed data
D <- length(prawn.data$ty)

# Generate posterior predictive draws using the 'post_pred_sn' function
predicted.counts.imperfect.case <- post_pred_sn(
  n = 13,                 # Number of animals per group
  D = D,                  # Number of groups
  m = 5,                  # Number of bins for discretization (if applicable)
  alpha = unlist(MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$alpha_sample),
  beta = unlist(MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$beta_sample),
  sensitivity = unlist(MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$sensitivity_sample)
)

# Extract simulated and observed count data
ty_rep_imperfect_case <- predicted.counts.imperfect.case$ty_rep
ty <- prawn.data$ty

# Check proportions of observed vs simulated counts (for model diagnostics)
prop.table(table(ty))                      # Observed data distribution
prop.table(table(ty_rep_imperfect_case))   # Simulated data distribution

# ----------------------------------------------------------------------------------#
# Posterior Predictive Check Visualization
# ----------------------------------------------------------------------------------#

# Generate barplots comparing observed and replicated counts
# (Option 1: Package version)
ppc.barplot.imperfect.case <- ppc_barplot(
  ty = ty,
  yrep = ty_rep_imperfect_case,
  trial = 13
)

# (Option 2: Custom version with fixed Y-axis scale)
ppc.barplot.imperfect.case <- ppc_barplot_fixed_y_scale(
  ty = ty,
  yrep = ty_rep_imperfect_case,
  trial = 13
)

# Access components of the posterior predictive plot object
ppc.barplot.imperfect.case$p.counts            # Raw counts plot
ppc.barplot.imperfect.case$p.prop              # Proportions plot
ppc.barplot.imperfect.case$p.prop.zero         # Zero-inflation plot
ppc.barplot.imperfect.case$p.zero.non.zero     # Combined plot
ppc.barplot.imperfect.case$prop.zero.non.zero  # Zero vs non-zero proportion values

# Customize plot appearance for the final figure
ppc.barplot.imperfect.case$p.zero.non.zero <-
  ppc.barplot.imperfect.case$p.zero.non.zero +
  ggtitle(NULL) +
  xlab(expression("Infected Groups (" * tilde(t)[yi] * ")")) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14, face = "plain"),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

# Save final plot to file
ggsave(
  filename = "PPC_Prawn_imperfect_Case_sn_0.80.png",
  plot = ppc.barplot.imperfect.case$p.zero.non.zero,
  width = 8,
  height = 4.5,
  units = "in"
)



# ----------------------------------------------------------------------------------
# Posterior Predictive Check: Imperfect Sensitivity with Uncertainty (Uniform(0.75, 0.85))
# Model: Exact Model using Metropolis-Hastings Algorithm
# Sensitivity Prior: Uniform(0.75, 0.85)
# Data: Imported Dataset
# Goal: Assess model fit using posterior predictive simulations under uncertain test sensitivity
# ----------------------------------------------------------------------------------

# Load posterior samples from MCMC output
load("MH.alpha.mu.sigma.20.unknown.sn.80.Prawn.Rdata")

# Extract number of test groups from observed data
D <- length(prawn.data$ty)

# Generate posterior predictive replicates under sensitivity uncertainty
predicted.counts.imperfect.uniform.case <- post_pred_sn(
  n = 13,   # Number of animals per group
  D = D,    # Number of groups
  m = 5,    # Number of bins (if relevant to the function)
  alpha = unlist(MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$alpha_sample),
  beta = unlist(MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$beta_sample),
  sensitivity = unlist(MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$sensitivity_sample)
)

# Extract simulated and observed group-level infection counts
ty_rep_imperfect_unif_case <- predicted.counts.imperfect.uniform.case$ty_rep
ty <- prawn.data$ty

# Summary: Proportion tables for comparison (optional diagnostics)
prop.table(table(ty))                        # Observed
prop.table(table(ty_rep_imperfect_unif_case)) # Simulated

# ----------------------------------------------------------------------------------#
# Posterior Predictive Visualization
# ----------------------------------------------------------------------------------#

# Generate PPC plots using a custom function with fixed Y-axis scaling
ppc.barplot.imperfect.uniform.case <- ppc_barplot_prawn_fixed_y_scale(
  ty = ty,
  yrep = ty_rep_imperfect_unif_case,
  trial = 13
)

# Access various plot components for reporting or diagnostics
ppc.barplot.imperfect.uniform.case$p.counts
ppc.barplot.imperfect.uniform.case$p.prop
ppc.barplot.imperfect.uniform.case$p.prop.zero
ppc.barplot.imperfect.uniform.case$p.zero.non.zero
ppc.barplot.imperfect.uniform.case$prop.zero.non.zero

# Customize final PPC plot appearance (zero vs non-zero infected groups)
ppc.barplot.imperfect.uniform.case$p.zero.non.zero <-
  ppc.barplot.imperfect.uniform.case$p.zero.non.zero +
  ggtitle(NULL) +
  xlab(expression("Infected Groups (" * tilde(t)[yi] * ")")) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14, face = "plain"),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black")
  )

# ----------------------------------------------------------------------------------#
# Save Final Plot
# ----------------------------------------------------------------------------------#

ggsave(
  filename = "PPC_Prawn_imperfect_Case_sn_Unif_0.80_revised.png",  # Output file name
  plot = ppc.barplot.imperfect.uniform.case$p.zero.non.zero,       # Plot to save
  width = 8,
  height = 4.5,
  units = "in",
  dpi = 300
)



