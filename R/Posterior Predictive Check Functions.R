#' Posterior Predictive Simulation for Grouped Beta-Binomial Models under perfect test
#'
#' Simulates posterior predictive responses \code{ty} under grouped beta-binomial models,
#' using either the approximate model (e.g., HMC-style with Beta(alpha, beta/m)) or
#' the exact model (e.g., MH-style with transformation via group pooling).
#'
#' @param n Integer. Number of trials per group (e.g., pool size).
#' @param alpha Numeric vector. Posterior samples of the alpha parameter (length equal to number of MCMC draws).
#' @param beta Numeric vector. Posterior samples of the beta parameter (same length as \code{alpha}).
#' @param D Integer. Number of distinct pools or groups (e.g., villages, domains).
#' @param m Integer. Number of individuals per group (used for transformation under exact model).
#' @param approximate.model Logical. If \code{TRUE}, use the approximate model (assumes prevalence ~ Beta(alpha, beta/m));
#'                          if \code{FALSE}, use the exact model with transformation \code{1 - (1 - p)^m}.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{\code{p_rep}}{A matrix of simulated prevalence probabilities (dimensions: draws × D).}
#'   \item{\code{ty_rep}}{A matrix of simulated binomial responses (dimensions: draws × D).}
#' }
#'
#' @examples
#' \dontrun{
#' alpha_samples <- rnorm(1000, 0.005, 0.2)
#' beta_samples <- rnorm(1000, 6, 0.5)
#' result <- post.pred.aprox.model(n = 13, alpha = alpha_samples, beta = beta_samples,
#'                                 D = 40, m = 5, approximate.model = TRUE)
#' }
#'
#' @export
post_pred_aprox_model <- function(n, alpha, beta, D, m, approximate.model = TRUE) {
  # alpha and beta are vectors
  # "approximate.model=TRUE" indicates whether the MCMC outputs are from Approximate model (HMC)
  # "approximate.model=FALSE" indicates whether the MCMC outputs are from Exact model (MH)

  no.sim <- length(alpha)  # Number of simulations (rows)

  # Initialize matrices
  p_matrix <- matrix(NA, nrow = no.sim, ncol = D)
  phi_i_matrix <- matrix(NA, nrow = no.sim, ncol = D)
  ty_matrix <- matrix(NA, nrow = no.sim, ncol = D)

  # Loop over each simulation
  for (i in 1:no.sim) {
    if (approximate.model) {
      p_matrix[i, ] <- rbeta(D, alpha[i], beta[i] / m)  # Batch-specific mean prevalence
      ty_matrix[i, ] <- rbinom(D, n, p_matrix[i, ])
    } else {
      p_matrix[i, ] <- rbeta(D, alpha[i], beta[i])  # Batch-specific mean prevalence
      phi_i_matrix[i, ] <- 1 - (1 - p_matrix[i, ])^m  # Adjusted probabilities
      ty_matrix[i, ] <- rbinom(D, n, phi_i_matrix[i, ])
    }
  }

  list(p_rep = p_matrix, ty_rep = ty_matrix)
}


#' Posterior Predictive Simulation for Grouped Beta-Binomial Model with Imperfect Sensitivity
#'
#' Simulates posterior predictive responses under a grouped beta-binomial model
#' with known test sensitivity applied to group-level prevalence estimates.
#' The model includes a transformation to account for pooling and imperfect test sensitivity.
#'
#' @param n Integer. Number of trials per group (e.g., pool size).
#' @param alpha Numeric vector. Posterior samples of the alpha parameter (length = number of MCMC draws).
#' @param beta Numeric vector. Posterior samples of the beta parameter (same length as \code{alpha}).
#' @param D Integer. Number of groups or domains (e.g., villages, pools).
#' @param m Integer. Number of individuals per pool (used in the transformation of prevalence).
#' @param sensitivity Numeric vector. Posterior samples of test sensitivity (same length as \code{alpha}).
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{p_rep}}{Matrix of simulated latent prevalences (dimensions: draws × D).}
#'   \item{\code{phi_i_matrix}}{Matrix of adjusted group-level detection probabilities
#'                              after accounting for pooling and sensitivity.}
#'   \item{\code{ty_rep}}{Matrix of simulated observed binomial counts (dimensions: draws × D).}
#' }
#'
#' @details
#' For each MCMC draw, prevalence in each group is simulated from a Beta distribution,
#' transformed to account for pooling as \code{1 - (1 - p)^m}, and then scaled by
#' the test sensitivity. The observed count is simulated using a Binomial distribution
#' with this adjusted probability.
#'
#' @examples
#' \dontrun{
#' alpha_samples <- rnorm(1000, 2, 0.1)
#' beta_samples <- rnorm(1000, 4, 0.2)
#' sn_samples <- runif(1000, 0.75, 0.80)
#' result <- post.pred.sn(n = 13, alpha = alpha_samples, beta = beta_samples,
#'                        D = 40, m = 5, sensitivity = sn_samples)
#' }
#'
#' @export
post_pred_sn <- function(n, alpha, beta, D, m, sensitivity) {
  # alpha, beta and sensitivity are vectors
  no.sim <- length(alpha)  # Number of simulations (rows)

  # Initialize matrices
  p_matrix <- matrix(NA, nrow = no.sim, ncol = D)
  phi_i_matrix <- matrix(NA, nrow = no.sim, ncol = D)
  ty_matrix <- matrix(NA, nrow = no.sim, ncol = D)

  # Loop over each simulation
  for (i in 1:no.sim) {
    # Simulate individual-level prevalence
    p_matrix[i, ] <- rbeta(D, alpha[i], beta[i])

    # Adjust for pooling and imperfect sensitivity
    phi_i_matrix[i, ] <- sensitivity[i] * (1 - (1 - p_matrix[i, ])^m)

    # Simulate observed test results
    ty_matrix[i, ] <- rbinom(D, n, phi_i_matrix[i, ])
  }

  list(p_rep = p_matrix, phi_i_matrix = phi_i_matrix, ty_rep = ty_matrix)
}



#' Posterior Predictive Check: Barplot for Observed vs. Simulated Counts and Proportions
#'
#' Generates diagnostic barplots for observed vs. posterior predictive simulated counts
#' and proportions in grouped binomial models, tailored for considered import case study.
#' Includes plots with 95% credible intervals, observed values, and side-by-side
#' comparisons of zero and non-zero group outcomes.
#'
#' @param ty Integer vector. Observed counts (e.g., number of infected bags per batch).
#' @param yrep Matrix. Simulated count outcomes from the posterior predictive distribution.
#'             Each row corresponds to one posterior draw, and columns to group-level outcomes.
#' @param trial Integer. Maximum number of possible positives (i.e., number of trials or test replicates per batch).
#'
#' @return A named list of ggplot objects and summary data frames:
#' \describe{
#'   \item{\code{count.summary}}{Summary statistics (mean, median, 95% CI) for count-based PPC.}
#'   \item{\code{prop.summary}}{Summary statistics for proportion-based PPC.}
#'   \item{\code{p.counts}}{Side-by-side barplots for counts (zero vs. non-zero) with 95% CI.}
#'   \item{\code{p.prop}}{Proportional version of \code{p.counts}.}
#'   \item{\code{p.prop.zero}}{Barplots comparing zero vs. non-zero proportions with CIs.}
#'   \item{\code{count.zero.summary}}{Summary of zero and non-zero count frequencies.}
#'   \item{\code{p.zero.non.zero}}{Dual-axis plot showing zero and non-zero count comparisons.}
#'   \item{\code{prop.zero.non.zero}}{Dual-axis plot showing zero and non-zero proportions.}
#' }
#'
#' @details
#' The function is designed for posterior predictive checking in models with grouped data,
#' such as pooled testing studies. It provides separate visualizations for:
#' - Raw count distributions (zero vs. non-zero infected groups)
#' - Proportions of groups by infection count
#' - Combined dual-axis plots showing comparisons between zero and non-zero groups
#'
#' Zero counts are handled separately to highlight sparsity in grouped outcomes.
#'
#' @examples
#' \dontrun{
#' # Simulate observed data and posterior draws
#' ty <- c(0, 1, 0, 2, 0, 3, 1)
#' yrep <- matrix(sample(0:3, 700, replace = TRUE), nrow = 100)
#' trial <- 3
#' plots <- ppc.barplot.prawn(ty, yrep, trial)
#'
#' # View the main count comparison plot
#' print(plots$p.counts)
#' }
#'
#' @import ggplot2 ggpubr
#' @export
ppc_barplot <- function(ty,yrep,trial){

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
    geom_col(aes(y = y2.observed*1/ratio.zero.nonzero), fill = "red", alpha = 0.5, width = 0.4) +  # Second y variable (y2), scaled and bars
    # Error bars for y1.observed (zero counts)
    geom_errorbar(aes(ymin = y1.025, ymax = y1.975), width = 0.2, color = "blue") +
    # Error bars for y2.observed (non-zero counts) - scaled
    geom_errorbar(aes(ymin = y2.025 * 1/ratio.zero.nonzero, ymax = y2.975 * 1/ratio.zero.nonzero), width = 0.2, color = "red") +

    scale_y_continuous(
      name = "Number of zero counts",         # Primary y-axis
      sec.axis = sec_axis(~.*ratio.zero.nonzero, name = "Number of non-zero counts")  # Secondary y-axis
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

  list(count.summary=count.summary,prop.summary=prop.summary,p.counts=p.counts,
       p.prop=p.prop,p.prop.zero=p.prop.zero,count.zero.summary=count.zero.summary,
       p.zero.non.zero=p.zero.non.zero,prop.zero.non.zero=prop.zero.non.zero)
}
