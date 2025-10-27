# R Functions

#' Generate Random Samples from a Truncated Beta Distribution
#'
#' Draws random samples from a Beta distribution and truncates them at a specified cutoff value.
#' If a sampled value is less than the cutoff, it is replaced with 0.
#'
#' @param n Integer. Number of samples to draw.
#' @param shape1 Numeric. First shape parameter of the Beta distribution.
#' @param shape2 Numeric. Second shape parameter of the Beta distribution.
#' @param cutoff Numeric between 0 and 1. Values sampled below this threshold are set to 0.
#'
#' @return A numeric vector of length \code{n}, containing truncated Beta samples.
#'
#' @examples
#' set.seed(123)
#' rtrbeta(10, shape1 = 2, shape2 = 5, cutoff = 0.2)
#'
#' @export
rtrbeta <- function(n, shape1, shape2, cutoff) {

  # Step 1: Draw contamination prevalence from Beta distribution
  beta_samples <- rbeta(n, shape1, shape2)

  # Step 2: Use cutoff for assigning design prevalence
  samples <- ifelse(beta_samples >= cutoff, beta_samples, 0)

  return(samples)
}

#' Compute the Logit of a Probability
#'
#' This function computes the logit (log-odds) transformation of a probability value.
#'
#' @param p A numeric value or vector strictly between 0 and 1. Represents the probability.
#'
#' @return A numeric value (or vector) representing the logit of the input probability.
#' @examples
#' logit(0.7)
#' logit(c(0.3, 0.5, 0.7))
#'
#' @export
logit <- function(p) {
  if (any(p <= 0 | p >= 1)) {
    stop("All p values must be between 0 and 1 (exclusive).")
  }
  return(log(p / (1 - p)))
}


#' Compute the Inverse Logit (Logistic Function)
#'
#' This function computes the inverse logit transformation, also known as the logistic function.
#' It maps real-valued inputs to the (0, 1) interval.
#'
#' @param x A numeric value or vector. Typically a log-odds or any real number.
#'
#' @return A numeric value (or vector) representing the inverse logit of the input.
#' @examples
#' inverse_logit(0.8472979)  # Approximately 0.7
#' inverse_logit(c(-1, 0, 1))
#'
#' @export
inverse_logit <- function(x) {
  exp_x <- exp(x)
  return(exp_x / (1 + exp_x))
}


#' Metropolis-Hastings Sampler with Truncated Beta Distribution
#'
#' Runs a Metropolis–Hastings (MH) MCMC algorithm to sample from a posterior
#' distribution involving a truncated Beta distribution. The function supports
#' several parameterizations, and optionally models test sensitivity and specificity
#' with prior distributions.
#'
#' @param par Initial parameter values (numeric vector).
#' @param cutoff Truncation threshold for the Beta distribution (numeric).
#'   cutoff = 0 leads to the simple beta-binomial model.
#' @param initial Starting values for the parameters (numeric vector).
#' @param alpha Logical. Whether to include `alpha` in the model.
#' @param beta Logical. Whether to include `beta` in the model.
#' @param mu Logical. Whether to include `mu` in the model.
#' @param rho Logical. Whether to include `rho` in the model.
#' @param G Integer. Number of quantiles from the beta distribution.
#' @param R Integer. Number of MCMC iterations.
#' @param sigma Proposal standard deviation for each parameter (numeric vector).
#' @param num_chains Number of MCMC chains to run (default is 3).
#' @param burnin Number of initial iterations to discard (default is `R * 0.5`).
#' @param thin Thinning interval for MCMC samples (default is 5).
#' @param seed A vector of integers for random number seeds (one per chain).
#' @param trail Integer. Number of trials (e.g., number of group-tests conducted per batch).
#' @param ty2 Integer vector of observed successes (e.g., number of groups tested positive).
#' @param wt Numeric vector of weights for each observation.
#' @param group.size Integer. Number of individuals tested per group (pool size).
#' @param sensitivity Logical. Whether to consider test sensitivity in the model.
#' @param specificity Logical. Whether to consider test specificity in the model.
#' @param sensitivity.range Optional numeric vector of length 2 indicating [min, max] for sensitivity
#'   (e.g., sensitivity ranges over U(0.65, 0.75)).
#' @param specificity.range Optional numeric vector of length 2 indicating [min, max] for specificity
#'   (e.g., specificity ranges over U(0.90, 0.99)).
#'
#' @importFrom MASS mvrnorm
#'
#' @return A list with four elements:
#' \describe{
#'   \item{target.parameters}{List of posterior samples for each parameter across chains.}
#'   \item{parameters}{List of MCMC samples for each chain.}
#'   \item{log.likelihood.values}{List of log-likelihood values for each iteration.}
#'   \item{log.likelihood.yi}{List of log-likelihood values for each observation at each iteration.}
#' }
#'
#' @details
#' Implements a flexible Bayesian MCMC sampler with support for multiple parameterizations
#' (via alpha, beta, mu, and rho). It includes likelihood computation using truncated
#' beta-distributed prevalence and allows modeling imperfect test sensitivity/specificity
#' using logit-transformed priors. The truncated beta-binomial model with cutoff k = 0
#' reduces to a standard beta-binomial model.
#'
#' @examples
#' \dontrun{
#' result <- MH_Sampler_BB(
#'   par = c(log(0.005), logit(7e-04)),
#'   cutoff = 0.2,
#'   initial = par,
#'   alpha = TRUE, beta = TRUE, mu = FALSE, rho = FALSE,
#'   G = 500, R = 1000,
#'   sigma = c(0.20, 0.20),
#'   num_chains = 4,
#'   burnin = 200, thin = 5,
#'   seed = c(123, 456, 789, 102),
#'   trail = 13,
#'   ty2 = 0:13,
#'   wt = c(2815,9,10,6,1,3,2,0,1,2,1,0,0,0),
#'   group.size = 5,
#'   sensitivity = FALSE,
#'   specificity = FALSE
#' )
#' }
#'
#' @export
MH_Sampler_BB <- function(par,cutoff,initial,alpha,beta,mu,rho,G,R,sigma,num_chains=3,burnin,thin=5,seed,trail,ty2,wt,group.size,sensitivity,specificity,
                                 sensitivity.range=NULL,specificity.range=NULL){

  # Provide seed number for each chain
  # cutoff: Truncation point
  target.parameters <- vector("list", num_chains)
  parameters <- vector("list", num_chains)
  log.likelihood.values <- vector("list", num_chains)
  log.likelihood.yi <- vector("list", num_chains)

  logit <- function(p) {
    if (p <= 0 || p >= 1) {
      stop("p must be between 0 and 1 (exclusive).")
    }
    return(log(p / (1 - p)))
  }

  inverse_logit <- function(x) {
    exp_x <- exp(x)
    return(exp_x / (1 + exp_x))
  }

  # Handle burnin
  if (is.null(burnin)) {
    burnin <- floor(R * 0.5)  # Default to 50% of R
  }

  if (burnin >= R) {
    stop("`burnin` must be less than the total number of iterations `R`.")
  }

  # Handle initial values
  if (missing(initial) || is.null(initial) || length(initial) == 0 || all(is.na(initial))) {
    initial <- par
  }


  prior.log.density.improper.uniform <- function(u){0} # Used for mu, alpha, beta and rho
  prior.log.density.logit <- function(u){u-2*log(1+exp(u))} # Used for sensitivity and specificity



  b= trail; Nbar=group.size; d=length(ty2)

  length.par <- length(par)

  for (chain in 1:num_chains) {


    par.matrix <- matrix(NA, R + 1, length.par)

    par.matrix[1, ] <- initial+ rnorm(length.par, 0, 1) # starting values
    set.seed(seed[chain])

    # Store Some other results

    loglik.value.results <- matrix(NA, R + 1, 3) # LL.old, LL.new, acceptance

    loglik.value.new.observation <- matrix(NA, R + 1, d) # For each y_i

    for (r in 1:R) {

      proposal <- par.matrix[r, ] + rnorm(length.par, 0, sigma)

      if (length(sigma) == 1) {
        proposal <- par.matrix[r, ] + rnorm(length.par, 0, sigma) # For scalar sigma
      } else if (is.vector(sigma) && length(sigma) == length.par) {
        proposal <- par.matrix[r, ] + mvrnorm(1, mu = rep(0, length.par), Sigma = diag(sigma^2)) # For vectore sigma
      }
      else if (is.matrix(sigma) && nrow(sigma) == length.par && ncol(sigma) == length.par) {
        # full covariance matrix
        proposal <- par.matrix[r, ] + mvrnorm(1, mu = rep(0, length.par), Sigma = sigma)
      } else {
        stop("sigma must be scalar, vector of length = length(par), or covariance matrix.")
      }


      if(alpha==TRUE & beta==TRUE) {
        # Compute shape parameters
        shape1_old <- exp(par.matrix[r, 1])
        shape2_old <- exp(par.matrix[r, 2])
        shape1_new <- exp(proposal[1])
        shape2_new <- exp(proposal[2])
      }

      if(alpha==TRUE & mu==TRUE) {
        # Compute shape parameters
        shape1_old <- exp(par.matrix[r, 1])
        shape2_old <- exp(par.matrix[r, 1]) * (1 - inverse_logit(par.matrix[r, 2])) / inverse_logit(par.matrix[r, 2])
        shape1_new <- exp(proposal[1])
        shape2_new <- exp(proposal[1]) * (1 - inverse_logit(proposal[2])) / inverse_logit(proposal[2])
      }

      if(mu==TRUE & rho==TRUE) {
        mu.old <- inverse_logit(par.matrix[r, 1])
        mu.new <- inverse_logit(proposal[1])

        rho.old <- inverse_logit(par.matrix[r, 2])
        rho.new <- inverse_logit(proposal[2])

        # Compute shape parameters
        shape1_old <- (1-rho.old)/rho.old*mu.old
        shape1_new <- (1-rho.new)/rho.new*mu.new
        shape2_old <- (1-rho.old)/rho.old*(1-mu.old)
        shape2_new <- (1-rho.new)/rho.new*(1-mu.new)
      }

      # Ensure shape parameters are positive
      shape1_old <- pmax(1e-10, shape1_old)
      shape2_old <- pmax(1e-10, shape2_old)
      shape1_new <- pmax(1e-10, shape1_new)
      shape2_new <- pmax(1e-10, shape2_new)


      # Compute quantiles
      p.i.vals.old <- qbeta(c(1:G) / (G + 1), shape1_old, shape2_old)
      p.i.vals.new <- qbeta(c(1:G) / (G + 1), shape1_new, shape2_new)

      p.i.vals.old <- p.i.vals.old * (p.i.vals.old>=cutoff)
      p.i.vals.new <- p.i.vals.new * (p.i.vals.new>=cutoff)

      phi.i.vals.old <- pmax(0, (1 - (1 - p.i.vals.old)^Nbar))
      phi.i.vals.new <- pmax(0, (1 - (1 - p.i.vals.new)^Nbar))

      if (sensitivity==TRUE & specificity==FALSE){
        sensitivity.min <- sensitivity.range[1]
        sensitivity.SD <- sensitivity.range[2]-sensitivity.range[1]

        # Calculation of sensitivity
        sensitivity_old <- sensitivity.min+sensitivity.SD*inverse_logit(par.matrix[r,length.par]) #  Delta=function(u)=inverse_logit(u)
        sensitivity_new <- sensitivity.min+sensitivity.SD*inverse_logit(proposal[length.par])

        # Compute phi values: With sensitivity
        phi.i.vals.old <- pmax(0, sensitivity_old*phi.i.vals.old)
        phi.i.vals.new <- pmax(0, sensitivity_new*phi.i.vals.new)

      }

      if (sensitivity==TRUE & specificity==TRUE){


        # Calculation of sensitivity
        sensitivity.min <- sensitivity.range[1]
        sensitivity.SD <- sensitivity.range[2]-sensitivity.range[1]

        sensitivity_old <- sensitivity.min+sensitivity.SD*inverse_logit(par.matrix[r,(length.par-1)]) #  Delta=function(u)=inverse_logit(u)
        sensitivity_new <- sensitivity.min+sensitivity.SD*inverse_logit(proposal[(length.par-1)])

        # Calculation of Specificity

        specificity.min <- specificity.range[1]
        specificity.SD <- specificity.range[2]-specificity.range[1]

        specificity_old <- specificity.min+specificity.SD*inverse_logit(par.matrix[r,(length.par)]) #  Delta=function(u)=inverse_logit(u)
        specificity_new <- specificity.min+specificity.SD*inverse_logit(proposal[(length.par)])

        # Compute phi values: With sensitivity and specificity
        phi.i.vals.old <- pmax(0, sensitivity_old*phi.i.vals.old+(1-specificity_old)*(1-phi.i.vals.old))
        phi.i.vals.new <- pmax(0, sensitivity_new*phi.i.vals.new+(1-specificity_new)*(1-phi.i.vals.new))

      }

      # Initialize log-likelihoods
      loglik.old <- 0
      loglik.new <- 0

      for (i in 1:d) {
        loglik.value.new.observation[r,i] <- wt[i] * log(mean(dbinom(x = ty2[i], size = b, prob = phi.i.vals.new)) + 1e-10)
        loglik.old <- loglik.old + wt[i] * log(mean(dbinom(x = ty2[i], size = b, prob = phi.i.vals.old)) + 1e-10)
        loglik.new <- loglik.new + wt[i] * log(mean(dbinom(x = ty2[i], size = b, prob = phi.i.vals.new)) + 1e-10)
      }


      if (sensitivity==FALSE & specificity==FALSE){

        prior.log.density.old <- prior.log.density.improper.uniform(par.matrix[r,1])+prior.log.density.improper.uniform(par.matrix[r,2])
        prior.log.density.new <- prior.log.density.improper.uniform(proposal[1])+prior.log.density.improper.uniform(proposal[2])

      }


      if (sensitivity==TRUE & specificity==FALSE){
        prior.log.density.old <- prior.log.density.improper.uniform(par.matrix[r,1])+prior.log.density.improper.uniform(par.matrix[r,2])+prior.log.density.logit(par.matrix[r,length.par])
        prior.log.density.new <- prior.log.density.improper.uniform(proposal[1])+prior.log.density.improper.uniform(proposal[2])+prior.log.density.logit(proposal[length.par])

      }


      if (sensitivity==TRUE & specificity==TRUE){
        prior.log.density.old <- prior.log.density.improper.uniform(par.matrix[r,1])+prior.log.density.improper.uniform(par.matrix[r,2])+prior.log.density.logit(par.matrix[r,(length.par-1)])+prior.log.density.logit(par.matrix[r,(length.par)])
        prior.log.density.new <- prior.log.density.improper.uniform(proposal[1])+prior.log.density.improper.uniform(proposal[2])+prior.log.density.logit(proposal[(length.par-1)])+prior.log.density.logit(proposal[(length.par)])
      }



      # Metropolis-Hastings acceptance criterion


      loglik.value.new <- loglik.new + prior.log.density.new
      loglik.value.old <- loglik.old + prior.log.density.old

      #    if (is.finite(loglik.new) && is.finite(loglik.old) && runif(1) <= exp(loglik.new + prior.log.density.new - loglik.old - prior.log.density.old)) {
      if (is.finite(loglik.new) && is.finite(loglik.old) && runif(1) <= exp(loglik.value.new - loglik.value.old)) {
        par.matrix[r + 1, ] <- proposal
        accepted <- TRUE
      } else {
        par.matrix[r + 1, ] <- par.matrix[r, ]
        accepted <- FALSE
      }

      loglik.value.results[r,] <- c(loglik.value.old,loglik.value.new,accepted)

    } # Loop Finish

    # Listing results

    par.matrix <- par.matrix[c((burnin+1):R),]
    sampled.par <- seq(1,dim(par.matrix)[1],by=thin)
    par.matrix <- par.matrix[sampled.par,]

    loglik.value.results <- loglik.value.results[c((burnin+1):R),]
    loglik.value.results <- loglik.value.results[sampled.par,]

    loglik.value.new.observation <- loglik.value.new.observation[c((burnin+1):R),]
    loglik.value.new.observation <- loglik.value.new.observation[sampled.par,]


    if (alpha==TRUE & beta==TRUE){

      alpha_sample <- exp(par.matrix[,1])
      beta_sample <- exp(par.matrix[,2])

      mu_sample <- alpha_sample/(alpha_sample+beta_sample)

      out <-  list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample)
    }

    if (alpha==TRUE & beta==TRUE & sensitivity==TRUE){
      alpha_sample <- exp(par.matrix[,1])
      beta_sample <- exp(par.matrix[,2])
      sensitivity_sample <- sensitivity.min+sensitivity.SD*inverse_logit(par.matrix[,3])
      mu_sample <- alpha_sample/(alpha_sample+beta_sample)
      out <-  list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample,sensitivity_sample=sensitivity_sample)
    }



    if (alpha==TRUE & mu==TRUE){
      alpha_sample <- exp(par.matrix[,1])
      mu_sample <- inverse_logit(par.matrix[,2])

      beta_sample <- alpha_sample*(1-mu_sample)/mu_sample
      out <- list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample)
    }


    if (alpha==TRUE & mu==TRUE & sensitivity==TRUE){
      alpha_sample <- exp(par.matrix[,1])
      mu_sample <- inverse_logit(par.matrix[,2])
      sensitivity_sample <- sensitivity.min+sensitivity.SD*inverse_logit(par.matrix[,3])
      beta_sample <- alpha_sample*(1-mu_sample)/mu_sample
      out <- list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample,sensitivity_sample=sensitivity_sample)
    }

    if (mu==TRUE & rho==TRUE){
      mu_sample <- inverse_logit(par.matrix[,1])
      rho_sample <- inverse_logit(par.matrix[,2])

      alpha_sample <- (1-rho_sample)/rho_sample*mu_sample
      beta_sample <- (1-rho_sample)/rho_sample*(1-mu_sample)
      out <-  list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample)
    }


    if (mu==TRUE & rho==TRUE & sensitivity==TRUE){
      mu_sample <- inverse_logit(par.matrix[,1])
      rho_sample <- inverse_logit(par.matrix[,2])

      alpha_sample <- (1-rho_sample)/rho_sample*mu_sample
      beta_sample <- (1-rho_sample)/rho_sample*(1-mu_sample)

      sensitivity_sample <- sensitivity.min+sensitivity.SD*inverse_logit(par.matrix[,3])

      out <- list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample,sensitivity_sample=sensitivity_sample)
    }




    if (alpha==TRUE & beta==TRUE & sensitivity==TRUE & specificity==TRUE){

      alpha_sample <- exp(par.matrix[,1])
      beta_sample <- exp(par.matrix[,2])
      mu_sample <- alpha_sample/(alpha_sample+beta_sample)

      sensitivity_sample <- sensitivity.min+sensitivity.SD*inverse_logit(par.matrix[,3])
      specificity_sample <- specificity.min+specificity.SD*inverse_logit(par.matrix[,4])

      out <- list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample,sensitivity_sample=sensitivity_sample,specificity_sample=specificity_sample)
    }


    if (alpha==TRUE & mu==TRUE & sensitivity==TRUE & specificity==TRUE){

      alpha_sample <- exp(par.matrix[,1])
      mu_sample <- inverse_logit(par.matrix[,2])
      beta_sample <- alpha_sample*(1-mu_sample)/mu_sample

      sensitivity_sample <- sensitivity.min+sensitivity.SD*inverse_logit(par.matrix[,3])
      specificity_sample <- specificity.min+specificity.SD*inverse_logit(par.matrix[,4])
      out <-  list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample,sensitivity_sample=sensitivity_sample,specificity_sample=specificity_sample)
    }


    if (mu==TRUE & rho==TRUE & sensitivity==TRUE & specificity==TRUE){

      mu_sample <- inverse_logit(par.matrix[,1])
      rho_sample <- inverse_logit(par.matrix[,2])

      alpha_sample <- (1-rho_sample)/rho_sample*mu_sample
      beta_sample <- (1-rho_sample)/rho_sample*(1-mu_sample)

      sensitivity_sample <- sensitivity.min+sensitivity.SD*inverse_logit(par.matrix[,3])
      specificity_sample <- specificity.min+specificity.SD*inverse_logit(par.matrix[,4])

      out <- list(alpha_sample=alpha_sample,beta_sample=beta_sample,mu_sample=mu_sample,sensitivity_sample=sensitivity_sample,specificity_sample=specificity_sample)

    }  # Complete task for one chain

    # Store the chain
    target.parameters[[chain]] <- out
    parameters[[chain]] <- par.matrix
    log.likelihood.values[[chain]] <- loglik.value.results
    log.likelihood.yi[[chain]] <- loglik.value.new.observation
  }

  # Get the names of all parameters in the first chain (assuming all chains have the same structure)
  param_names <- names(target.parameters[[1]])

  # Initialize an empty list to hold combined chains for each parameter
  combined_list <- list()

  # Loop through each parameter and combine its values across all chains
  for (param in param_names) {
    # Combine the values from all chains for the current parameter
    combined_list[[param]] <- lapply(target.parameters, function(chain) chain[[param]])
  }


  list(target.parameters=combined_list, parameters=parameters,log.likelihood.values=log.likelihood.values,log.likelihood.yi=log.likelihood.yi)
}


#' Estimate Expected Leakage and Probability of Leakage
#'
#' Computes the expected leakage and probability of leakage under a Beta-Binomial model.
#' Supports vector/scalar inputs or posterior samples (lists of chains).
#' Optionally applies a truncation cutoff on the Beta distribution tail.
#'
#' @param alpha Numeric vector or list of posterior chains for alpha (Beta distribution shape1).
#' @param beta Numeric vector or list of posterior chains for beta (Beta distribution shape2).
#' @param B Numeric scalar. Upper bound on leakage parameter.
#' @param b Numeric scalar. Lower bound on leakage parameter.
#' @param M Numeric scalar. Total number of items/entities.
#' @param m Numeric scalar. Number of items/entities considered in partial leakage.
#' @param sensitivity Numeric scalar, vector, or list of sensitivities. If list, treated as posterior samples. If sensitivity=1, exact formula will be used for probability of leakage.
#' @param cutoff Numeric scalar in [0, 1) or NULL. If specified, truncates calculations to [cutoff, 1] interval.
#'
#' @return A list containing inputs and computed quantities:
#' \describe{
#'   \item{B}{Upper leakage bound (input).}
#'   \item{b}{Lower leakage bound (input).}
#'   \item{cutoff}{Cutoff threshold for truncation, or NULL.}
#'   \item{M}{Total number of items (input).}
#'   \item{m}{Number of items considered in partial leakage (input).}
#'   \item{alpha}{Input alpha parameter(s).}
#'   \item{beta}{Input beta parameter(s).}
#'   \item{mu}{Mean of Beta distribution(s), alpha / (alpha + beta).}
#'   \item{sensitivity}{Input sensitivity parameter(s).}
#'   \item{E.L}{Expected leakage estimate(s). If posterior samples, a nested list; otherwise numeric vector.}
#'   \item{P.L}{Probability of leakage estimate(s). If posterior samples, a nested list; otherwise numeric vector.}
#' }
#'
#' @details
#' The function estimates leakage under a Beta-Binomial model by integrating
#' over truncated or full Beta distributions, depending on whether a cutoff is specified.
#' It supports posterior chains by accepting lists of alpha, beta, and sensitivity values.
#' The expected leakage is calculated as a
#' ratio of incomplete Beta integrals, weighted by group sizes.
#' Two formulations for probability of leakage are used depending
#' on whether `sensitivity` equals 1 (exact) or differs from 1 (approximation).

#' @examples
#' # Scalar/vector input example
#' estimate_leakage(alpha = 0.005, beta = 4.5, B = 8000, b = 13, M = 40, m = 5, sensitivity = 0.8)
#'
#' # Posterior sample input example (lists of numeric vectors)
#' alpha_chain <- list(c(0.004, 0.005, 0.006), c(0.0035, 0.0042, 0.0055))
#' beta_chain <- list(c(3, 5, 6), c(4, 3, 5))
#' sensitivity_chain <- list(c(0.75, 0.8, 0.85), c(0.80, 0.75, 0.75))
#' estimate_leakage(alpha_chain, beta_chain, B = 8000, b = 13, M = 40, m = 5, sensitivity_chain,cutoff=0.01)
#' estimate_leakage(alpha_chain, beta_chain, B = 8000, b = 13, M = 40, m = 5, sensitivity = 0.8,cutoff=0.01)
#'
#' @export
estimate_leakage <- function(alpha, beta, B, b, M, m, sensitivity, cutoff = NULL) {


  approximate <- if (is.list(sensitivity)) {
    # For posterior chains, assume approximation if any value < 1
    any(sapply(sensitivity, function(s) any(s < 1)))
  } else {
    any(sensitivity < 1, na.rm = TRUE)
  }


  # ---- Internal helper: Incomplete Beta over [cutoff, 1] ---- #
  incomplete_lbeta_k1 <- function(alpha, beta, k, min_log = -1e5) {
    if (any(k < 0 | k >= 1)) stop("cutoff k must be in [0,1)")
    if (any(alpha <= 0 | beta <= 0)) stop("alpha, beta must be > 0")

    log_B_complete <- lbeta(alpha, beta)
    tail_prob <- 1 - pbeta(k, alpha, beta)
    log_tail_prob <- ifelse(tail_prob > 0, log(tail_prob), min_log)

    log_B_complete + log_tail_prob
  }

  # ---- Core: Expected & Probability of Leakage (Truncated) ---- #
  compute_truncated <- function(a, b_, s) {
    bmDelta <- b * m * s

    # Expected Leakage
    num_trunc <- exp(incomplete_lbeta_k1(a + 1, b_ + bmDelta, k = cutoff))
    denom_trunc <- exp(incomplete_lbeta_k1(a, b_, k = 0))
    beta_ratio <- num_trunc / denom_trunc
    E_L <- (B - b) * M * beta_ratio

    # Probability of Leakage
    beta_mDelta <- b_ / (m * s)

    if (approximate) {
      min_phi_m <- s * (1 - (1 - cutoff)^m)
      num1 <- exp(incomplete_lbeta_k1(a, beta_mDelta + b, k = min_phi_m))
      denom1 <- exp(incomplete_lbeta_k1(a, beta_mDelta, k = 0))
      term1 <- num1 / denom1
    } else {
      num1 <- exp(incomplete_lbeta_k1(a, b_ + m * b, k = cutoff))
      denom1 <- exp(incomplete_lbeta_k1(a, b_, k = 0))
      term1 <- num1 / denom1
    }

    num2 <- exp(incomplete_lbeta_k1(a, b_ + M * B, k = cutoff))
    denom2 <- exp(incomplete_lbeta_k1(a, b_, k = 0))
    term2 <- num2 / denom2

    P_L <- term1 - term2
    return(list(E_L = E_L, P_L = P_L))
  }

  # ---- Core: Expected & Probability of Leakage (Standard, non-truncated) ---- #

  compute_standard <- function(a, b_, s) {
    if (m == M) {
      # Expected Leakage
      E_L <- (B - b) * M * exp(lbeta(a + 1, b_ + b * M * s) - lbeta(a, b_))

      # Probability of Leakage
      if (approximate) {
        P_L <- exp(lbeta(a, b_ / (M * s) + b) - lbeta(a, b_ / (M * s))) -
          exp(lbeta(a, b_ + M * B) - lbeta(a, b_))
      } else {
        P_L <- exp(lbeta(a, b_ + M*b) - lbeta(a, b_ )) -
          exp(lbeta(a, b_ + M * B) - lbeta(a, b_))
      }

    } else {
      # Expected Leakage
      E_L <- (B - b) * M * exp(lbeta(a + 1, b_ + b * m * s) - lbeta(a, b_))

      # Probability of Leakage
      if (approximate) {
        P_L <- exp(lbeta(a, b_ / (m * s) + b) - lbeta(a, b_ / (m * s))) -
          exp(lbeta(a, b_ + M * B) - lbeta(a, b_))
      } else {
        P_L <- exp(lbeta(a, b_ +m * b) - lbeta(a, b_ )) -
          exp(lbeta(a, b_ + M * B) - lbeta(a, b_))
      }
    }

    return(list(E_L = E_L, P_L = P_L))
  }




  # ---- Compute mu ---- #
  mu <- if (is.list(alpha)) {
    Map(function(a, b_) a / (a + b_), alpha, beta)
  } else {
    alpha / (alpha + beta)
  }

  # ---- Main logic: handle lists or scalars ---- #
  if (is.list(alpha)) {
    # posterior samples case
    if (is.list(sensitivity)) {
      # posterior samples for alpha, beta, and sensitivity
      results <- Map(function(a_chain, b_chain, s_chain) {
        Map(function(a, b_, s) {
          if (!is.null(cutoff)) {
            compute_truncated(a, b_, s)
          } else {
            compute_standard(a, b_, s)
          }
        }, a_chain, b_chain, s_chain)
      }, alpha, beta, sensitivity)
    } else {
      # posterior samples for alpha and beta, scalar sensitivity
      results <- Map(function(a_chain, b_chain) {
        Map(function(a, b_) {
          if (!is.null(cutoff)) {
            compute_truncated(a, b_, sensitivity)
          } else {
            compute_standard(a, b_, sensitivity)
          }
        }, a_chain, b_chain)
      }, alpha, beta)
    }
    # return(results)

    # Extract E.L and P.L from results
    # results is a nested list: list of chains of lists with E_L and P_L
    # E_L <- lapply(results, function(chain) lapply(chain, `[[`, "E_L"))
    # P_L <- lapply(results, function(chain) lapply(chain, `[[`, "P_L"))


    E_L <- lapply(results, function(chain) {
      sapply(chain, `[[`, "E_L")  # Returns a numeric vector
    })

    P_L <- lapply(results, function(chain) {
      sapply(chain, `[[`, "P_L")  # Returns a numeric vector
    })


    return(list(
      B = B,
      b = b,
      cutoff = cutoff,
      M = M,
      m = m,
      alpha = alpha,
      beta = beta,
      mu = mu,
      sensitivity = sensitivity,
      E.L = E_L,
      P.L = P_L
    ))
  } else {
    # vector/scalar alpha, beta, sensitivity

    # Ensure sensitivity is same length or scalar
    if (length(sensitivity) == 1) {
      sensitivity_vec <- rep(sensitivity, length(alpha))
    } else {
      sensitivity_vec <- sensitivity
    }

    out <- mapply(function(a, b_, s) {
      if (!is.null(cutoff)) {
        compute_truncated(a, b_, s)
      } else {
        compute_standard(a, b_, s)
      }
    }, alpha, beta, sensitivity_vec, SIMPLIFY = FALSE)

    E_L <- unlist(lapply(out, `[[`, "E_L"))
    P_L <- unlist(lapply(out, `[[`, "P_L"))

    return(list(
      B = B,
      b = b,
      cutoff = cutoff,
      M = M,
      m = m,
      alpha = alpha,
      beta = beta,
      mu = mu,
      sensitivity = sensitivity,
      E.L = E_L,
      P.L = P_L
    ))
  }

}

#' Summarize MCMC Chains with Diagnostics
#'
#' This function provides a summary of MCMC output, including mean, standard deviation,
#' quantiles, Gelman-Rubin diagnostics (R-hat), effective sample size (n_eff), and
#' efficiency ratios for multiple chain inputs. It supports both single-chain and multiple-chain (chained) inputs.
#'
#' @param ... One or more MCMC outputs. For multiple chains, each parameter should be
#'   passed as a list of numeric vectors (one per chain). For a single chain, pass each
#'   parameter as a single numeric vector.
#' @param varnames A character vector of parameter names (must match the number of parameters).
#'
#' @return A list containing:
#' \describe{
#'   \item{summary_mcmc}{A data frame with summary statistics for each parameter:
#'   Mean, SD, quantiles (2.5%, 25%, 50%, 75%, 97.5%), R-hat, n_eff, and n_eff_ratio.}
#'   \item{n_eff_chain}{A list of effective sample sizes per chain.}
#'   \item{n_eff_ratio_chain}{A list of effective sample size ratios per chain.}
#' }
#'
#' @examples
#' # --- Example 1: Single-chain input (vectors for alpha and beta) ---
#' set.seed(123)
#' alpha <- rnorm(1000, mean = 0.005, sd = 0.1)
#' beta <- rnorm(1000, mean = 5, sd = 0.2)
#'
#' result_single <- summary_mcmc(alpha, beta, varnames = c("alpha", "beta"))
#' print(result_single$summary_mcmc)
#'
#'
#' # --- Example 2: Multiple chains (list of vectors per parameter) ---
#' set.seed(456)
#' n_iter <- 1000
#'
#' # Simulate three chains for each parameter
#' alpha_chains <- list(
#'   rnorm(n_iter, mean = 0.005, sd = 0.1),
#'   rnorm(n_iter, mean = 0.005, sd = 0.1),
#'   rnorm(n_iter, mean = 0.005, sd = 0.1)
#' )
#'
#' beta_chains <- list(
#'   rnorm(n_iter, mean = 5, sd = 0.5),
#'   rnorm(n_iter, mean = 5, sd = 0.5),
#'   rnorm(n_iter, mean = 5, sd = 0.5)
#' )
#'
#' # Run the summary
#' result_multi <- summary_mcmc(alpha_chains, beta_chains, varnames = c("alpha", "beta"))
#' print(result_multi$summary_mcmc)
#'
#' @importFrom coda mcmc mcmc.list gelman.diag effectiveSize
#' @export
summary_mcmc <- function(..., varnames) {
  library(coda)

  # Capture arguments and detect structure
  raw_inputs <- list(...)

  is_chained <- is.list(raw_inputs[[1]]) && all(sapply(raw_inputs[[1]], is.numeric))

  if (is_chained) {
    ##### ----- MULTIPLE CHAINS ----- #####

    vector_lists <- raw_inputs
    no.chain <- length(vector_lists[[1]])

    combined_chains <- vector("list", length = no.chain)
    for (i in 1:no.chain) {
      combined_chains[[i]] <- do.call(cbind, lapply(vector_lists, `[[`, i))
    }

    chain_sizes <- sapply(combined_chains, nrow)
    mcmc_chains_parameters <- mcmc.list(lapply(combined_chains, mcmc))

    for (i in seq_along(mcmc_chains_parameters)) {
      varnames(mcmc_chains_parameters[[i]]) <- varnames
    }

    total_param_lengths_chain <- sum(chain_sizes)

    # Gelman-Rubin diagnostics
    R_hat_transformed_parameter <- gelman.diag(mcmc_chains_parameters, confidence = 0.95)
    R_hat <- c(R_hat_transformed_parameter$psrf[, 1])

    # Effective sample size
    n_eff <- effectiveSize(mcmc_chains_parameters)
    n_eff_chain <- lapply(mcmc_chains_parameters, effectiveSize)
    n_eff_ratio <- n_eff / total_param_lengths_chain
    n_eff_ratio_chain <- lapply(n_eff_chain, function(eff) eff / chain_sizes[1])  # assumes same length

    # Summary
    summ_param <- summary(mcmc_chains_parameters)
    summary_mcmc <- cbind(
      summ_param$statistics[, c("Mean", "SD")],
      summ_param$quantiles,
      R_hat = R_hat,
      n_eff = n_eff,
      n_eff_ratio = n_eff_ratio
    )

    rownames(summary_mcmc) <- varnames

    return(list(
      summary_mcmc = summary_mcmc,
      n_eff_chain = n_eff_chain,
      n_eff_ratio_chain = n_eff_ratio_chain
    ))

  } else {
    ##### ----- SINGLE CHAIN (VECTOR) ----- #####

    matrix_params <- do.call(cbind, raw_inputs)
    colnames(matrix_params) <- varnames

    mcmc_single <- as.mcmc(matrix_params)

    # Use internal summary method to avoid conflicts
    summ_param <- coda:::summary.mcmc(mcmc_single)

    n_eff <- NA
    n_eff_ratio <- NA

    summary_mcmc <- cbind(
      summ_param$statistics[, c("Mean", "SD")],
      summ_param$quantiles,
      R_hat = NA,
      n_eff = n_eff,
      n_eff_ratio = n_eff_ratio
    )

    rownames(summary_mcmc) <- varnames

    return(list(
      summary_mcmc = summary_mcmc,
      n_eff_chain = list(n_eff),
      n_eff_ratio_chain = list(n_eff_ratio)
    ))
  }
}


#' Create MCMC Array for Bayesplot Visualization
#'
#' Converts a list of MCMC chains (typically from a Metropolis-Hastings algorithm)
#' into a 3D array suitable for plotting using the `bayesplot` package.
#'
#' @param mcmc.list.chain A list of 4 matrices or data frames, each representing
#'   one MCMC chain. Each matrix should have parameters as columns and iterations as rows.
#' @param varnames A character vector of variable names corresponding to the columns
#'   in each chain (e.g., \code{c("log(alpha)", "logit(mu)")}). These will be used
#'   as parameter names in the output array.
#'
#' @return A 3-dimensional array with dimensions (iterations × chains × parameters),
#'   suitable for use with `bayesplot` functions such as \code{mcmc_trace()} or
#'   \code{mcmc_dens_overlay()}.
#'
#' @importFrom coda mcmc mcmc.list
#' @examples
#' # Example usage (assuming 4 chains of posterior samples):
#' \dontrun{
#'   chains_list <- list(chain1, chain2, chain3, chain4)
#'   param_names <- c("log(alpha)", "logit(mu)")
#'   mcmc_array <- create_mcmc_BP(chains_list, param_names)
#'   bayesplot::mcmc_trace(mcmc_array)
#' }
#'
#' @export
create_mcmc_BP <- function(mcmc.list.chain, varnames) {
  library(coda)

  # Convert to 'mcmc.list' object
  mcmc_chains <- mcmc.list(
    mcmc(mcmc.list.chain[[1]]),
    mcmc(mcmc.list.chain[[2]]),
    mcmc(mcmc.list.chain[[3]]),
    mcmc(mcmc.list.chain[[4]])
  )

  # Assign variable names
  varnames(mcmc_chains) <- varnames

  # Convert to array
  mcmc_array <- as.array(mcmc_chains)

  # Ensure proper shape: (iterations x chains x parameters)
  if (length(dim(mcmc_array)) > 2) {
    mcmc_array <- aperm(mcmc_array, c(1, 3, 2))  # reorder dims
  } else {
    mcmc_array <- array(mcmc_array, dim = c(dim(mcmc_array), 1))
  }

  # Define dimension names
  dimnames(mcmc_array) <- list(
    iterations = seq_len(dim(mcmc_array)[1]),
    chains     = paste0("chain", seq_len(dim(mcmc_array)[2])),
    parameters = varnames
  )

  return(mcmc_array)
}




