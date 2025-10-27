#' Log-Likelihood for Grouped Beta-Binomial Model with Test Error
#'
#' Computes the (log-)likelihood for the grouped Beta-Binomial (BB) model, accounting for test sensitivity and specificity,
#' and optionally computes expected and observed leakage under perfect specificity.
#'
#' @param par A numeric vector of log-parameters to estimate. The interpretation depends on which parameters are fixed or missing:
#'   - Typically includes `log(alpha)`, `log(beta)`, `log(sensitivity)`, and `log(specificity)` as needed.
#' @param ty A numeric vector of grouped test outcomes (e.g., number of positives in each group).
#' @param freq A numeric vector of frequencies (how many times each group outcome occurred).
#' @param theta Truncation parameter. Use `Inf` for the standard truncated Beta distribution.
#' @param mu Optional. Mean parameter of the Beta distribution. If supplied, `beta` is derived as `(1 - mu) / mu * alpha`.
#' @param alpha Optional. Alpha shape parameter of the Beta distribution.
#' @param sensitivity Optional. Test sensitivity. If missing, it is estimated via `par`.
#' @param specificity Optional. Test specificity. If missing, it is estimated via `par`.
#' @param b Integer. Number of individuals tested in each group (binomial trials).
#' @param R Integer. Number of grid points used to approximate the truncated integral (default 1000).
#' @param m Integer. Number of subsampled individuals per pool.
#' @param M Integer. Number of individuals per super-pool (used for leakage calculation).
#' @param deviance Logical. If `TRUE`, returns deviance (i.e., `-2 * loglikelihood`); otherwise returns log-likelihood.
#' @param subtract Numeric. A constant to subtract from the likelihood (useful for nested model comparison).
#' @param B Integer. Total number of positive individuals in super-pooling (for leakage calculation).
#' @param mushape Logical. If `TRUE`, uses a parameterization where `mu` and `shape` are estimated instead of `alpha` and `beta`.
#' @param alpha.mu Logical. If `TRUE`, uses a parameterization with `alpha` and `mu` as parameters.
#' @param leakage Logical. If `TRUE`, returns a named vector including leakage probability (`P.L`) and expected leakage (`E.L`).
#'
#' @return A numeric value representing the log-likelihood (or deviance) of the model. If `leakage = TRUE`, a named numeric vector is returned:
#'   - `deviance`: The deviance or log-likelihood.
#'   - `P.L`: Probability of leakage.
#'   - `E.L`: Expected number of leaked positives.
#'
#' @examples
#' Example usages of `loglik_group_bb` with different parameter missingness cases.
#' These examples illustrate how to fit the Beta-Binomial model under different
#' assumptions about which parameters are missing and need to be estimated.
#'
#' @examples
#' # 1. alpha and mu missing: sensitivity and specificity given
#' MLE.BB.m1 <- optim(c(0,0), loglik_group_bb, ty=ty, freq=freq, b=13,
#'                 theta=Inf, m=5, M=40, deviance=FALSE,
#'                 control=list(reltol=1e-12, fnscale=-1),
#'                 R=1e4, hessian=FALSE, sensitivity=1, specificity=1)
#' MLE.BB.m1$mu <- exp(MLE.BB.m1$par[1]) / sum(exp(MLE.BB.m1$par))
#' MLE.BB.m1$alpha <- exp(MLE.BB.m1$par[1])
#' MLE.BB.m1$beta <- exp(MLE.BB.m1$par[2])
#'
#' # 2. alpha missing (mu provided): sensitivity and specificity given
#' MLE.BB.m2 <- optim(c(0.005), loglik_group_bb, ty=ty, freq=freq, b=13,
#'                       theta=Inf, mu=0.0007, sensitivity=1, specificity=1,
#'                       m=5, M=40, deviance=FALSE,
#'                       control=list(factr=1e-12, fnscale=-1, trace=TRUE),
#'                       R=1e4, hessian=FALSE, method="L-BFGS-B",
#'                       lower=c(-Inf), upper=c(Inf))
#' MLE.BB.m2$alpha <- exp(MLE.BB.m2$par[1])
#' MLE.BB.m2$beta <- MLE.BB.m2$alpha * (1 - 0.0007) / 0.0007

#' # 3. mu missing (alpha provided): sensitivity and specificity given
#' MLE.BB.m3 <- optim(c(0.0007), loglik_group_bb, ty=ty, freq=freq, b=13,
#'                       theta=Inf, alpha=0.005, sensitivity=1, specificity=1,
#'                       m=5, M=40, deviance=FALSE,
#'                       control=list(factr=1e-12, fnscale=-1, trace=TRUE),
#'                       R=1e4, hessian=FALSE, method="L-BFGS-B",
#'                       lower=c(-Inf), upper=c(Inf))
#' MLE.BB.m3$beta <- exp(MLE.BB.m3$par[1])
#'
#' # 4. alpha, mu, and sensitivity missing: specificity given
#' MLE.BB.m4 <- optim(c(0,0,0), loglik_group_bb, ty=ty, freq=freq, b=13,
#'                       theta=Inf, m=5, M=40, specificity=1, deviance=FALSE,
#'                       control=list(factr=1e-12, fnscale=-1, trace=TRUE),
#'                       R=1e4, hessian=FALSE, method="L-BFGS-B",
#'                       lower=c(log(0.0001), log(0.0001), log(0.50)),
#'                       upper=c(Inf, Inf, 0))
#' MLE.BB.m4$mu <- exp(MLE.BB.m4$par[1]) / sum(exp(MLE.BB.m4$par))
#' MLE.BB.m4$alpha <- exp(MLE.BB.m4$par[1])
#' MLE.BB.m4$beta <- exp(MLE.BB.m4$par[2])
#' MLE.BB.m4$sensitivity <- exp(MLE.BB.m4$par[3])
#'
#' # 5. alpha, mu, and specificity missing: : sensitivity given
#' MLE.BB.m5 <- optim(c(log(0.005), log(3.0), log(0.5)), loglik_group_bb,
#'                       ty=ty, freq=freq, b=13, theta=Inf, m=5, M=40,
#'                       sensitivity=1, deviance=FALSE,
#'                       control=list(factr=1e-12, fnscale=-1, trace=TRUE),
#'                       R=1e4, hessian=FALSE, method="L-BFGS-B",
#'                       lower=c(log(0.0001), log(0.0001), log(0.99)),
#'                       upper=c(Inf, Inf, 0))
#' MLE.BB.m5$mu <- exp(MLE.BB.m5$par[1]) / sum(exp(MLE.BB.m5$par))
#' MLE.BB.m5$alpha <- exp(MLE.BB.m5$par[1])
#' MLE.BB.m5$beta <- exp(MLE.BB.m5$par[2])
#' MLE.BB.m5$specificity <- exp(MLE.BB.m5$par[3])
#'
#' # 6. alpha, mu, sensitivity, and specificity are all missing
#' MLE.BB.m6 <- optim(c(0,0,0,0), loglik_group_bb,
#'                       ty=ty, freq=freq, b=13, theta=Inf, m=5, M=40,
#'                       deviance=FALSE,
#'                       control=list(factr=1e-12, fnscale=-1, trace=TRUE),
#'                       R=1e4, hessian=FALSE, method="L-BFGS-B",
#'                       lower=c(log(0.0001), log(0.0001), log(0.5), log(0.99)),
#'                       upper=c(Inf, Inf, 0, 0))
#' MLE.BB.m6$mu <- exp(MLE.BB.m6$par[1]) / sum(exp(MLE.BB.m6$par))
#' MLE.BB.m6$alpha <- exp(MLE.BB.m6$par[1])
#' MLE.BB.m6$beta <- exp(MLE.BB.m6$par[2])
#' MLE.BB.m6$sensitivity <- exp(MLE.BB.m6$par[3])
#' MLE.BB.m6$specificity <- exp(MLE.BB.m6$par[4])

#' @export

loglik_group_bb <- function(par,ty,freq,theta,mu,alpha,beta,sensitivity, specificity,b,R=1000,m,M,deviance=FALSE,
                            subtract=0,B,mushape=FALSE,alpha.mu=FALSE,leakage=FALSE){

  keep <- freq != 0
  ty <- ty[keep]
  freq <- freq[keep]

  Nbar=m

  # === Determine alpha and beta depending on inputs ===
  missing.alpha <- missing(alpha)
  missing.beta  <- missing(beta)
  missing.mu    <- missing(mu)

  if (!missing.mu) {
    # If mu is provided
    if (missing.alpha) {
      alpha <- exp(par[1])
    }
    beta <- (1 - mu) / mu * alpha
  } else if (!missing.beta && missing.alpha) {
    # If beta is provided but alpha is missing → estimate alpha
    alpha <- exp(par[1])
  } else if (missing.mu && missing.alpha && missing.beta) {
    # If both alpha and beta missing → estimate both
    alpha <- exp(par[1])
    beta  <- exp(par[2])
  } else if (!missing.alpha && missing.beta) {
    # If alpha is provided and beta is missing → estimate beta
    beta <- exp(par[1])
  }


  # For given mu: Missing alpha, sensitivity and specificity
  if(missing.alpha & !missing(mu) & missing(sensitivity) & missing(specificity)) {
    sensitivity <- exp(par[2])
    specificity <- exp(par[3])
  }

  # For given alpha: Missing mu, sensitivity and specificity
  if(!missing.alpha & missing(mu) & missing(sensitivity) & missing(specificity) ) {
    sensitivity <- exp(par[2])
    specificity <- exp(par[3])
  }

  if(missing.alpha & missing(mu) & !missing(sensitivity) & missing(specificity)) specificity <- exp(par[3])
  if(missing.alpha & missing(mu) & missing(sensitivity)  & !missing(specificity)) sensitivity <- exp(par[3])

  if(missing(sensitivity) & missing(specificity) & missing.alpha & missing(mu)) {
    sensitivity <- exp(par[3])
    specificity <- exp(par[4]) }



  if(mushape){
    mu <- 1-1/(1+exp(par[1]))
    shape <- exp(par[2])
    alpha <- shape*mu
    beta <- shape*(1-mu)
  }
  if(alpha.mu){
    alpha <- exp(par[1])
    mu <- 1 / (1+exp(par[2]))
    beta <- alpha*(1-mu)/mu
  }




  if(theta>0){
    pvals <- qbeta(c(1:R)/(R+1),alpha,beta)
    log.pvals <- log(pvals)
    if(theta==Inf){
      log.one.minus.phi.vals.perfect <- Nbar*Rmpfr::log1mexp(-log.pvals) # this is phi when no test errors
      log.phi.vals.perfect <- Rmpfr::log1mexp(-log.one.minus.phi.vals.perfect)
      one.minus.phi.vals.perfect <- exp(log.one.minus.phi.vals.perfect)
      phi.vals.perfect <- exp(log.phi.vals.perfect)

      phi.vals <- (1-specificity)*one.minus.phi.vals.perfect + sensitivity*phi.vals.perfect
      log.one.minus.phi.vals <- log(1-phi.vals)
      log.phi.vals <- log(phi.vals)
    }

    out <- 0
    for(j in c(1:length(ty))){
      if(ty[j]!=0) log.terms <- ty[j]*log.phi.vals +
          (b-ty[j])*log.one.minus.phi.vals
      if(ty[j]==0) log.terms <- b*log.one.minus.phi.vals
      K <- max(log.terms)-10
      out <- out + freq[j]*(log(mean(exp(log.terms-K)))+K)
    }
  }
  if(leakage){
    if(theta==Inf){
      # P.L <- exp(lbeta(alpha,beta/(Nbar*sensitivity)+b)-lbeta(alpha,beta/(Nbar*sensitivity)))-
      #   exp(lbeta(alpha,beta/(M*sensitivity)+B)-lbeta(alpha,beta/(M*sensitivity)))
      # P.L <- exp(lbeta(alpha,beta/(Nbar*sensitivity)+b)-lbeta(alpha,beta/(Nbar*sensitivity)))-
      #  exp(lbeta(alpha,beta/(M)+B)-lbeta(alpha,beta/(M)))
      P.L <- exp(lbeta(alpha, beta / (Nbar * sensitivity) + b) - lbeta(alpha, beta / (Nbar * sensitivity))) -
        exp(lbeta(alpha, beta + M*B) - lbeta(alpha, beta))
      E.L <- (B-b)*M*exp(lbeta(alpha+1,beta+b*Nbar*sensitivity)-lbeta(alpha,beta))
    }
  }
  if((!leakage)&deviance) return(-2*out-subtract)
  if((!leakage)&(!deviance)) return(out-subtract)
  if(leakage) return(c(deviance=-2*out,P.L=P.L,E.L=E.L))
}


#' Log-Likelihood Function for Group Testing with Test Error under minimum prevalence
#'
#' Computes the log-likelihood under a truncated Beta-Binomial model for grouped testing outcomes, with support for imperfect diagnostic tests.
#' This function supports flexible parameterizations (e.g., in terms of alpha, beta, mu) and includes extensions for truncated prevalence modeling,
#' numerical integration over latent prevalence values, and optional leakage estimation under imperfect sensitivity.
#' Suitable for use with optimization routines (e.g., \code{optim}) to estimate model parameters from grouped test outcome data.
#'
#' @param par A numeric vector of parameters (on log or logit scale) depending on what is fixed or estimated.
#' @param ty Observed number of positive group tests per replicate.
#' @param freq Frequency of each ty value.
#' @param theta Clustering parameter. Use `Inf` for considering random groups.
#' @param mu Optional; mean of the Beta prior for contamination probability.
#' @param alpha Optional; shape parameter of the Beta prior.
#' @param beta Optional; shape parameter of the Beta prior.
#' @param cutoff Truncation cutoff for the Beta distribution. `cutoff=0` indicates standard BB model.
#' @param sensitivity Test sensitivity.
#' @param specificity Test specificity.
#' @param b Number of groups per replicate.
#' @param R Number of integration points for numerical approximation.
#' @param M Group size.
#' @param m Pool-size used for group testing. This can be equal to M if all the group members are used for testing.
#' @param deviance Logical. If TRUE, returns deviance (-2 * log-likelihood).
#' @param subtract A constant to subtract from the log-likelihood.
#' @param B Reference number of groups per batch / replicate for leakage modeling.
#' @param mushape Logical. If TRUE, parameterize using mu and shape.
#' @param alpha.mu Logical. If TRUE, parameterize using alpha and mu.
#' @param leakage Logical. If TRUE, compute expected leakage quantities assuming imperfect sensitivity.
#'
#' @return Depending on arguments:
#' \itemize{
#'   \item If `leakage = FALSE` and `deviance = TRUE`, returns deviance.
#'   \item If `leakage = FALSE` and `deviance = FALSE`, returns log-likelihood.
#'   \item If `leakage = TRUE`, returns a named vector with deviance, expected leakage, and probability leakage.
#' }
#'
#' @details This function integrates over a truncated Beta distribution to compute the expected likelihood,
#' adjusted for misclassification due to imperfect test sensitivity and specificity. It supports various
#' parameterizations for alpha, beta, and mu, and incorporates numerical integration over latent prevalence.
#' Parameter configuration depends on which inputs are missing:
#' \itemize{
#' \item If \code{theta}, \code{mu}, and \code{alpha} are all missing, then \code{par = c(log(alpha), log(beta), log(theta))}.
#' \item If \code{theta} is missing, it is inferred as \code{exp(last element of par)}.
#' \item If \code{alpha} is missing, it is inferred as \code{exp(first element of par)}.
#' \item If \code{alpha} is provided, \code{beta} is taken as \code{exp(first element of par)} (when \code{mu} is missing).
#' \item If \code{mu} is supplied, then \code{beta = (1 - mu) / mu * alpha}.
#' }
#'
#' Users can directly call this function within optimization routines such as \code{optim()} to obtain maximum likelihood estimates.
#' The user should specify fixed parameters (e.g., mu or alpha) and pass remaining parameters through \code{par}.
#' For maximum flexibility, the function is not wrapped in a specific optimizer, allowing users to tailor optimization (e.g., convergence criteria, bounds, or transformation) based on their needs.
#' The output of \code{optim()} provides the parameter estimates on the log or logit scale, depending on the model specification.
#' Integration is performed over the latent prevalence distribution:
#' \itemize{
#' \item \code{R} controls the number of quantile points used for numeric integration over \code{p}.
#' }
#'
#' Expected parameter formats:
#' \itemize{
#' \item For given \code{alpha}: \code{par = c(beta, sensitivity, specificity)}
#' \item For given \code{mu}: \code{par = c(alpha, sensitivity, specificity)}
#' \item For given \code{beta}: \code{par = c(alpha, sensitivity, specificity)}
#' }
#'
#' @note
#' This implementation assumes \code{theta = Inf}, which indicates groups are formed randomly (uncluster).
#'
#' @section Leakage:
#' Leakage quantities (\code{E.L}, \code{P.L}) are calculated under the assumption of perfect specificity.
#' @importFrom Rmpfr log1mexp
#' @examples
#'freq = c(2815,9,10,6,1,3,2,0,1,2,1,0,0,0)
#'ty = c(0:13)
#'# For estimating log(alpha) and log(beta) for a cutoff=0, which reflects standard BB model
#'mle.est.cutoff.0 <- optim(c(0, 0), loglik_group_trbb,ty = ty,freq = freq,b=13,m=5,R = 1e4,theta = Inf,cutoff = 0,control = list(reltol = 1e-12, fnscale = -1))
#'exp(mle.est.cutoff.0$par)
#'# For estimating log(alpha) and log(beta) for a given cutt-off of 1%, which reflects truncated BB model
#'mle.est.cutoff.01 <- optim(c(0, 0), loglik_group_trbb,ty = ty,freq = freq,b=13,m=5,theta = Inf,R = 1e4,sensitivity = 1, specificity = 1,cutoff = 0.01,deviance = FALSE,control = list(reltol = 1e-12, fnscale = -1),hessian = FALSE)
#'exp(mle.est.cutoff.01$par)
#'# For estimating log(beta) for a given alpha and cutoff greater than 0, which reflects truncated BB model
#'beta.est.cutoff.01 <- optim(c(0),loglik_group_trbb,ty=ty,freq=freq,b=13,m=5,theta=Inf,R=1e4,alpha=0.005,cutoff = 0.01,sensitivity=1,specificity=1,
#'deviance=FALSE,method="Brent", control=list(reltol=1e-12,fnscale=-1), hessian=FALSE,lower = c(log(0.05)),upper = c(log(100)))
#'exp(beta.est.cutoff.01$par)
#'# For estimating log(alpha) for a given beta and cutoff greater than 0, which reflects truncated BB model
#'alpha.est.cutoff.01 <- optim(c(0),loglik_group_trbb,ty=ty,freq=freq,b=13,m=5,theta=Inf,R=1e4,beta=4.5,cutoff = 0.01,sensitivity=1,specificity=1,
#'deviance=FALSE,method="Brent", control=list(reltol=1e-12,fnscale=-1), hessian=FALSE,lower = c(log(0.0005)),upper = c(log(0.5)))
#'exp(alpha.est.cutoff.01$par)
#'# For estimating log(alpha) for a given mu and cutoff greater than 0, which reflects truncated BB model
#'# alpha missing but mu available
#'alpha.est.cutoff.01 <- optim(c(0),loglik_group_trbb,ty=ty,freq=freq,b=13,m=5,theta=Inf,R=1e4,mu=0.0005,cutoff = 0.01,sensitivity=1,specificity=1,
#'deviance=FALSE,method="Brent", control=list(reltol=1e-12,fnscale=-1), hessian=FALSE,lower = c(log(0.00005)),upper = c(log(0.05)))
#'exp(alpha.est.cutoff.01$par)
#'# For estimating leakage for a given alpha, beta, sensitivity and cutoff greater than 0, which reflects truncated BB model
#'  leakage.estimate <- loglik_group_trbb(c(log(0.005),log(4.5)),ty=ty,freq=freq,b=13,B=8000,m=5,M=40,cutoff=0.01,
#'  theta=Inf,leakage=TRUE,R=1e4,sensitivity=1.0,specificity=1.0)
#' leakage.estimate
#' @export
loglik_group_trbb <- function(par,ty,freq,theta,mu,alpha,beta,cutoff,sensitivity=NULL, specificity=NULL,b,R=1000,m,M,deviance=FALSE,
                                      subtract=0,B,mushape=FALSE,alpha.mu=FALSE,leakage=FALSE){

  keep <- freq != 0
  ty <- ty[keep]
  freq <- freq[keep]

  Nbar=m
  missing.alpha <- missing(alpha)

  if(missing(alpha)) alpha <- exp(par[1])
  if(missing(mu)&missing.alpha&missing(beta)) beta <- exp(par[2])
  if(missing(mu)&!missing.alpha) beta <- exp(par[1])
  if(!missing(beta) & missing(alpha)&missing(mu)) alpha <- exp(par[1])
  if(!missing(mu)&missing(beta)) beta <- (1-mu)/mu*alpha


  # At the start: ensure sensitivity and specificity have default values if NULL
  if (is.null(sensitivity)) sensitivity <- NA_real_
  if (is.null(specificity)) specificity <- NA_real_



  # For given mu: Missing alpha, beta, sensitivity and specificity
  if(missing.alpha & missing(beta) & !missing(mu) & missing(sensitivity) & missing(specificity)) {
    sensitivity <- exp(par[2])
    specificity <- exp(par[3])
  }

  # For given alpha: Missing beta, mu, sensitivity and specificity
  if(!missing.alpha & missing(beta) & missing(mu) & missing(sensitivity) & missing(specificity) ) {
    sensitivity <- exp(par[2])
    specificity <- exp(par[3])
  }

  # For given beta: Missing alpha, mu, sensitivity and specificity
  if(missing.alpha & !missing(beta) & missing(mu) & missing(sensitivity) & missing(specificity) ) {
    sensitivity <- exp(par[2])
    specificity <- exp(par[3])
  }

  if(missing.alpha & missing(beta) & missing(mu) & !missing(sensitivity) & missing(specificity)) specificity <- exp(par[3])
  if(missing.alpha & missing(beta) & missing(mu) & missing(sensitivity)  & !missing(specificity)) sensitivity <- exp(par[3])

#  if(missing(sensitivity) & missing(specificity) & missing.alpha & missing(mu) & missing(beta)) {
#    sensitivity <- exp(par[3])
#    specificity <- exp(par[4]) }

  # Case: infer sensitivity and specificity from par
  if (is.null(sensitivity) && is.null(specificity) && missing(alpha) && missing(mu) && missing(beta)) {
    sensitivity <- exp(par[3])
    specificity <- exp(par[4])
  }

  if(mushape){
    mu <- 1-1/(1+exp(par[1]))
    shape <- exp(par[2])
    alpha <- shape*mu
    beta <- shape*(1-mu)
  }
  if(alpha.mu){
    alpha <- exp(par[1])
    mu <- 1 / (1+exp(par[2]))
    beta <- alpha*(1-mu)/mu
  }

  # If after all this sensitivity or specificity are still NA, default to 1
  if (is.na(sensitivity)) sensitivity <- 1
  if (is.na(specificity)) specificity <- 1


  if(theta>0){


    cdf.cutoff <- pbeta(cutoff,alpha,beta) # Computes the cumulative probability up to the cutoff under the Beta(α, β) distribution.
    pvals <- qbeta(cdf.cutoff+(1-cdf.cutoff)*c(1:R)/(R+1), alpha, beta) #Generates R quantile points (like a grid) from the truncated Beta distribution over [cutoff,1].
    # pvals <- pvals * (pvals >= cutoff)
    log.pvals <- log(pvals)

    if(theta==Inf){
      log.one.minus.phi.vals.perfect <- Nbar*Rmpfr::log1mexp(-log.pvals) # this is phi when no test errors
      log.phi.vals.perfect <- Rmpfr::log1mexp(-log.one.minus.phi.vals.perfect)
      one.minus.phi.vals.perfect <- exp(log.one.minus.phi.vals.perfect)
      phi.vals.perfect <- exp(log.phi.vals.perfect)

      phi.vals <- (1-specificity)*one.minus.phi.vals.perfect + sensitivity*phi.vals.perfect
      log.one.minus.phi.vals <- log(1-phi.vals)
      log.phi.vals <- log(phi.vals)
    }


    out <- 0
    for(j in c(1:length(ty))){
      if(ty[j]!=0) log.terms <- ty[j]*log.phi.vals +
          (b-ty[j])*log.one.minus.phi.vals
      if(ty[j]==0) log.terms <- b*log.one.minus.phi.vals
      K <- max(log.terms)-10
      out <- out + freq[j]*(log((1-cdf.cutoff)*mean(exp(log.terms-K))+cdf.cutoff*(ty[j]==0)*exp(-K))+K)
      # # numerically stable approximation to the integral:𝐸_{𝜋 > cutoff}[P(Ty = ty | 𝜋)]

    }
  }

  if(leakage){
    if(theta==Inf){
      # P.L <- exp(lbeta(alpha,beta/(Nbar*sensitivity)+b)-lbeta(alpha,beta/(Nbar*sensitivity)))-
      #   exp(lbeta(alpha,beta/(M*sensitivity)+B)-lbeta(alpha,beta/(M*sensitivity)))
      # P.L <- exp(lbeta(alpha,beta/(Nbar*sensitivity)+b)-lbeta(alpha,beta/(Nbar*sensitivity)))-
      #   exp(lbeta(alpha,beta/(M)+B)-lbeta(alpha,beta/(M)))
      # E.L <- (B-b)*M*exp(lbeta(alpha+1,beta+b*Nbar*sensitivity)-lbeta(alpha,beta))
      E.L <- (B - b) * M * exp(lbeta(alpha + 1, beta + b * Nbar * sensitivity) - lbeta(alpha, beta))
      if(sensitivity==1){P.L <- exp(lbeta(alpha, beta + Nbar * b) - lbeta(alpha, beta)) -
        exp(lbeta(alpha, beta + B*M) - lbeta(alpha, beta))} else
      {P.L <- exp(lbeta(alpha, beta / (Nbar * sensitivity) + b) - lbeta(alpha, beta / (Nbar * sensitivity))) -
        exp(lbeta(alpha, beta + B*M) - lbeta(alpha, beta))}

    }
  }
  if((!leakage)&deviance) return(-2*out-subtract)
  if((!leakage)&(!deviance)) return(out-subtract)
  if(leakage) return(c(deviance=-2*out,P.L=P.L,E.L=E.L))
}







#' Fit Grouped Beta-Binomial (GroupedBB) Model to Grouped Data
#'
#' Fits a beta-binomial model to grouped count data using maximum likelihood estimation.
#' Estimates parameters alpha and beta unless provided, and supports leakage estimation
#' when `leakage = TRUE`.
#'
#' @param ty Numeric vector. The observed counts of "successes" in grouped data.
#' @param freq Numeric vector. Frequencies corresponding to each count in `ty`.
#' @param b Integer. Number of trials per group.
#' @param m Integer. Minimum group size or related parameter used in likelihood calculation.
#' @param M Integer, optional. Maximum group size or related parameter used in likelihood calculation. Required if `leakage = TRUE`.
#' @param B Integer, optional. Number of bootstrap or simulation iterations required if `leakage = TRUE`.
#' @param sensitivity Numeric, default 1. Sensitivity parameter for misclassification or detection.
#' @param specificity Numeric, default 1. Specificity parameter for misclassification or detection.
#' @param theta Numeric, default `Inf`. Overdispersion or shape parameter.
#' @param R Integer, default 10000. Number of repetitions or simulations used in likelihood approximation.
#' @param leakage Logical, default FALSE. Whether to include leakage in the model. Requires `M` and `B` if TRUE.
#' @param control List. Control parameters passed to `optim()`.
#' @param hessian Logical, default FALSE. Whether to return the Hessian matrix from optimization.
#'
#' @return A list containing estimated parameters `alpha`, `beta`, `mu` (mean parameter),
#'         sensitivity, specificity, log-likelihood value (`loglik`), and convergence status.
#'         If `hessian = TRUE`, also returns the Hessian matrix.
#'         If `leakage = TRUE`, includes leakage-related estimates `P.L` and `E.L`.
#'
#' @examples
#' freq <- c(2815, 9, 10, 6, 1, 3, 2, 0, 1, 2, 1, 0, 0, 0)
#' ty <- 0:13
#'
#' # Fit model without leakage
#' bb.m1 <- fit_GroupedBB(ty = ty, freq = freq, b = 13, m = 5, M = 40,
#'                        sensitivity = 1, specificity = 1)
#'
#' # Fit model with leakage (requires B)
#' bb.m2 <- fit_GroupedBB(ty = ty, freq = freq, b = 13, m = 5, M = 40, B = 8000,
#'                        sensitivity = 1, specificity = 1, leakage = TRUE)
#'
#' @export
fit_GroupedBB <- function(ty, freq, b, m,
                          sensitivity = 1, specificity = 1,
                          theta = Inf, R = 1e4,
                          control = list(reltol = 1e-12, fnscale = -1),
                          hessian = FALSE,
                          leakage = FALSE,
                          M = NULL, B = NULL) {

  if (leakage) {
    if (is.null(M)) stop("Parameter 'M' must be provided when leakage = TRUE")
    if (is.null(B)) stop("Parameter 'B' must be provided when leakage = TRUE")
  }

  par_init <- c(0, 0)  # initial guess for log(alpha), log(beta)

  fit <- optim(par_init, loglik_group_bb,
               ty = ty, freq = freq, b = b,
               m = m, M = M,
               theta = theta,
               sensitivity = sensitivity, specificity = specificity,
               control = control, R = R,
               hessian = hessian)

  alpha_hat <- exp(fit$par[1])
  beta_hat <- exp(fit$par[2])
  mu_hat <- alpha_hat / (alpha_hat + beta_hat)

  results <- list(alpha = alpha_hat,
                  beta = beta_hat,
                  mu = mu_hat,
                  sensitivity = sensitivity,
                  specificity = specificity,
                  loglik = fit$value,
                  convergence = fit$convergence)

  if (hessian) results$hessian <- fit$hessian

  if (leakage) {
    leakage_info <- loglik_group_bb(fit$par, ty = ty, freq = freq, b = b,
                                    theta = theta, m = m, M = M,
                                    B = B, deviance = FALSE, leakage = TRUE, R = R,
                                    sensitivity = sensitivity, specificity = specificity)

    results$P.L <- leakage_info["P.L"]
    results$E.L <- leakage_info["E.L"]
  }

  return(results)
}


#' Fit Truncated Grouped Beta-Binomial (trGroupedBB) Model to Grouped Data
#'
#' This function fits a truncated beta-binomial model to grouped count data using maximum likelihood estimation
#' for a given cutoff. It estimates parameters alpha and beta unless provided, and supports estimating leakage
#' when `leakage = TRUE`.
#'
#' @param ty Numeric vector. Observed counts.
#' @param freq Numeric vector. Frequencies for each count.
#' @param b Integer. Number of trials per group.
#' @param m Integer. Minimum group size or related parameter.
#' @param cutoff Numeric. Truncation cutoff value.
#' @param sensitivity Numeric, default 1. Sensitivity parameter.
#' @param specificity Numeric, default 1. Specificity parameter.
#' @param theta Numeric, default Inf. Overdispersion or shape parameter.
#' @param R Integer, default 1000. Number of quantiles for numerical integration.
#' @param control List. Control parameters passed to optim.
#' @param return_hessian Logical, default FALSE. Return Hessian matrix if TRUE.
#' @param leakage Logical, default FALSE. Whether to model leakage. Requires M and B if TRUE.
#' @param M Integer, optional. Required if leakage = TRUE.
#' @param B Integer, optional. Required if leakage = TRUE.
#'
#' @return A list containing estimated parameters `alpha`, `beta`, `mu` (mean parameter),
#'         sensitivity, specificity, log-likelihood value (`loglik`), and convergence status.
#'         If `hessian = TRUE`, also returns the Hessian matrix. If `leakage = TRUE`, includes
#'         leakage-related estimates `P.L` and `E.L`.
#'
#' @examples
#' freq <- c(2815, 9, 10, 6, 1, 3, 2, 0, 1, 2, 1, 0, 0, 0)
#' ty <- 0:13
#'
#' # Fit truncated grouped beta-binomial model with cutoff = 0 (no truncation)
#' trbb.m1 <- fit_trGroupedBB(ty, freq, b = 13, m = 5, cutoff = 0)
#' print(trbb.m1)
#'
#' # Fit truncated grouped beta-binomial model with cutoff = 0.01 and imperfect sensitivity
#' trbb.m2 <- fit_trGroupedBB(ty, freq, b = 13, m = 5, cutoff = 0.01,
#'                        sensitivity = 0.8, specificity = 1)
#' print(trbb.m2)
#' @export
fit_trGroupedBB <- function(ty, freq, b, m, cutoff,
                            sensitivity = 1, specificity = 1,
                            theta = Inf,
                            R = 1000,
                            control = list(reltol = 1e-12, fnscale = -1),
                            return_hessian = FALSE,
                            leakage = FALSE,
                            M = NULL,
                            B = NULL) {

  # Validate leakage-related parameters
  if (leakage) {
    if (is.null(M)) stop("Parameter 'M' must be provided when leakage = TRUE")
    if (is.null(B)) stop("Parameter 'B' must be provided when leakage = TRUE")
  }

  # Initial parameter guess for log(alpha) and log(beta)
  par_init <- c(0, 0)

  # Perform maximum likelihood estimation via optim
  fit <- optim(
    par = par_init,
    fn = loglik_group_trbb,
    ty = ty,
    freq = freq,
    b = b,
    m = m,
    cutoff = cutoff,
    sensitivity = sensitivity,
    specificity = specificity,
    theta = theta,
    R = R,
    deviance = FALSE,
    control = control,
    hessian = return_hessian
  )

  # Extract estimated parameters (on original scale)
  alpha_hat <- exp(fit$par[1])
  beta_hat  <- exp(fit$par[2])

  # Compile results
  results <- list(
    alpha = alpha_hat,
    beta = beta_hat,
    sensitivity = sensitivity,
    specificity = specificity,
    cutoff = cutoff,
    loglik = fit$value,
    convergence = fit$convergence
  )

  if (return_hessian) {
    results$hessian <- fit$hessian
  }

  # Compute leakage statistics if requested
  if (leakage) {
    leakage_info <- loglik_group_trbb(
      fit$par, ty = ty, freq = freq, b = b, m = m, cutoff = cutoff,
      theta = theta, M = M, B = B, deviance = FALSE, leakage = TRUE,
      R = R, sensitivity = sensitivity, specificity = specificity
    )
    results$P.L <- leakage_info["P.L"]
    results$E.L <- leakage_info["E.L"]
  }

  return(results)
}


#' Profile Likelihood Confidence Interval Estimation
#'
#' Computes a confidence interval for a parameter based on its profile log-likelihood.
#' The function identifies the maximum likelihood estimate (MLE) and determines the
#' interval bounds where the log-likelihood falls below the cutoff determined by a
#' chi-squared distribution with one degree of freedom.
#'
#' @param parameter Numeric vector of parameter values corresponding to the evaluated log-likelihoods.
#' @param logL Numeric vector of log-likelihood values for the corresponding parameters.
#' @param level Confidence level for the interval (default = 0.95).
#'
#' @details
#' The profile likelihood confidence interval is based on the likelihood ratio principle:
#' \deqn{2(\ell(\hat{\theta}) - \ell(\theta)) \sim \chi^2_1}
#' where \eqn{\ell(\hat{\theta})} is the maximum log-likelihood. The cutoff value for
#' the log-likelihood is computed as:
#' \deqn{\ell_{\text{cutoff}} = \ell(\hat{\theta}) - \frac{1}{2}\chi^2_{1, \text{level}}}
#'
#' The lower and upper bounds of the interval are identified as the parameter values
#' where the log-likelihood is closest to the cutoff, on either side of the MLE.
#'
#' @return A list containing:
#' \item{mle}{The maximum likelihood estimate (MLE).}
#' \item{mle.ll}{The lower confidence limit.}
#' \item{mle.ul}{The upper confidence limit.}
#' \item{level}{The specified confidence level.}
#' \item{value}{The log-likelihood cutoff value.}
#'
#' @examples
#' # Example: simple quadratic log-likelihood
#' param <- seq(0, 10, length.out = 100)
#' logL <- -0.5 * (param - 5)^2
#' pllfCIestimate(param, logL, level = 0.95)
#'
#' @export
pllfCIestimate <- function(parameter, logL, level = 0.95) {
  # Check input lengths
  if (length(parameter) != length(logL)) {
    stop("Parameter and log-likelihood vectors must be of equal length.")
  }

  # Remove NAs
  valid_idx <- which(!is.na(logL) & !is.na(parameter))
  parameter <- parameter[valid_idx]
  logL <- logL[valid_idx]

  # Find MLE and max log-likelihood
  max_logL <- max(logL)
  mle_idx <- which.max(logL)
  mle <- parameter[mle_idx]

  # Compute critical value from chi-squared
  chi_crit <- qchisq(level, df = 1)
  logL_cutoff <- max_logL - chi_crit / 2

  # Get CI bounds by minimizing distance from cutoff on each side of MLE
  lower_candidates <- which(parameter < mle)
  upper_candidates <- which(parameter > mle)

  # Handle edge cases
  mle.ll <- NA
  if (length(lower_candidates) > 0) {
    lower_diff <- abs(logL[lower_candidates] - logL_cutoff)
    mle.ll <- parameter[lower_candidates[which.min(lower_diff)]]
  }

  mle.ul <- NA
  if (length(upper_candidates) > 0) {
    upper_diff <- abs(logL[upper_candidates] - logL_cutoff)
    mle.ul <- parameter[upper_candidates[which.min(upper_diff)]]
  } else {
    mle.ul <- max(parameter)
  }

  # Return results
  list(
    mle = mle,
    mle.ll = mle.ll,
    mle.ul = mle.ul,
    level = level,
    value = logL_cutoff
  )
}
