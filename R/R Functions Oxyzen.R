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
#' Runs a Metropolis-Hastings (MH) MCMC algorithm to sample from a posterior distribution involving
#' a truncated Beta distribution. The function supports several parametrizations, and optionally models
#' test sensitivity and specificity with prior distributions.
#'
#' @param par Initial parameter values (numeric vector).
#' @param cutoff Truncation threshold for the Beta distribution (numeric).
#' @param initial Starting values for the parameters (numeric vector).
#' @param alpha Logical. Whether to include `alpha` in the model.
#' @param beta Logical. Whether to include `beta` in the model.
#' @param mu Logical. Whether to include `mu` in the model.
#' @param rho Logical. Whether to include `rho` in the model.
#' @param G Integer. Number of quantiles from beta distribution.
#' @param R Integer. Number of MCMC iterations.
#' @param sigma Proposal standard deviation for each parameter (numeric vector).
#' @param num_chains Number of MCMC chains to run (default is 3).
#' @param burnin Number of initial iterations to discard (default is `R * 0.5`).
#' @param thin Thinning interval for MCMC samples (default is 5).
#' @param seed A vector of integers for random number seed (one per chain).
#' @param trail Integer. Number of trials (say, number of group-test conducted per batch).
#' @param ty2 Integer. vector of observed successes (say, no. of groups tested positive).
#' @param wt Numeric vector of weights for each observation (say, no. of batches with no contamination).
#' @param group.size Integer. Number of individuals tested per group, in another way pool size.
#' @param sensitivity Logical. Whether to consider test sensitivity in model development.
#' @param specificity Logical. Whether to consider test specificity in model development.
#' @param sensitivity.range Optional numeric vector of length 2 indicating [min, max] for sensitivity. As for example, sensitivity ranges over U(0.65, 0.75)
#' @param specificity.range Optional numeric vector of length 2 indicating [min, max] for specificity. As for example, specificity. ranges over U(0.90, 0.99)
#'
#' @return A list containing:
#' \item{target.parameters}{List of posterior samples for each parameter across chains.}
#' \item{parameters}{List of MCMC samples for each chain.}
#' \item{log.likelihood.values}{List of log-likelihood values for each iteration.}
#' \item{log.likelihood.yi}{List of log-likelihood values for each observation at each iteration.}
#'
#' @details
#' The function implements a flexible Bayesian MCMC sampler with support for multiple parametrizations
#' (via alpha, beta, mu, and rho). It includes likelihood computation with binomial models using truncated
#' Beta-distributed prevalence and allows incorporating imperfect test sensitivity/specificity using
#' logit-transformed priors. The truncated beta-binomial model with cutoff $k=0$ is nothing but a standard beta-binomial model.
#'
#' @examples
#' \dontrun{
#' result <- MH.Sampler.Chain.TBB(
#'   par = c(0.1, 0.1),
#'   cutoff = 0.2,
#'   initial = c(0.1, 0.1),
#'   alpha = TRUE, beta = TRUE, mu = FALSE, rho = FALSE,
#'   G = 1000, R = 1000,
#'   sigma = c(0.20, 0.20),
#'   num_chains = 3,
#'   burnin = 5000, thin = 10,
#'   seed = c(123, 456, 789),
#'   trail = 13,
#'   ty2 = c(0:13),
#'   wt = c(2815,9,10,6,1,3,2,0,1,2,1,0,0,0),
#'   group.size = 5,
#'   sensitivity = FALSE,
#'   specificity = FALSE
#' )
#' }
#'
#' @export
MH.Sampler.Chain.TBB <- function(par,cutoff,initial,alpha,beta,mu,rho,G,R,sigma,num_chains=3,burnin=R*0.5,thin=5,seed,trail,ty2,wt,group.size,sensitivity,specificity,
                                 sensitivity.range=NULL,specificity.range=NULL){
  # Provide seed number for each chain
  # cutoff: Truncation point
  target.parameters <- vector("list", num_chains)
  parameters <- vector("list", num_chains)
  log.likelihood.values <- vector("list", num_chains)
  log.likelihood.yi <- vector("list", num_chains)


  prior.log.density.improper.uniform <- function(u){0} # Used for mu, alpha, beta and rho
  prior.log.density.logit <- function(u){u-2*log(1+exp(u))} # Used for sensitivity and specificity



  b= trail; Nbar=group.size;

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


#' Expected Leakage and Probability of leakage from Truncated Beta Distribution
#'
#' Computes the expected leakage and the probability of leakage from a  formulation using
#' a truncated Beta distribution. The function provides values from analytic forms and uses log-scale incomplete
#' Beta functions to avoid underflow. Analytic function are based on the perfect specificity. The probability of leakage is estimated using an approximate distribution for the group probability of being contaminated when sensitivity is imperfect. For perfect sensitivity, there is no approximate distribution required.
#'
#' @param alpha Numeric vector of shape parameters for the Beta distribution.
#' @param beta Numeric vector of shape parameters for the Beta distribution.
#' @param cutoff Numeric value in (0, 1). The truncation point for the Beta distribution.
#' @param B Integer. Total number of groups (say, bags) in a batch.
#' @param b Integer. Sampled number of groups (say, sample bags used in group-test) for a batch.
#' @param m Integer. Pool-size, i.e, number of sampled units used in group-test from each sampled group.
#' @param M Integer. Total number of units in a group (bag).
#' @param Delta Numeric scalar for sensitivity. Scale factor to adjust sensitivity (default = 1).
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{alpha}{Input vector of alpha values.}
#'   \item{beta}{Input vector of beta values.}
#'   \item{E_Li}{Expected leakage for each (alpha, beta) pair.}
#'   \item{P_Li}{Estimated probability of leakage for each (alpha, beta) pair.}
#'   \item{Note}{Estimation of expected and probability of leakage are under perfect specificty.}
#' }
#'
#' @details
#' This function evaluates integrals over the tail of the Beta distribution (i.e., from `cutoff` to 1)
#' using log-scale calculations to maintain numerical accuracy. The expected leakage is calculated as a
#' ratio of incomplete Beta integrals, weighted by group sizes. Two formulations for probability of leakage are used depending
#' on whether `Delta` equals 1 (exact) or differs from 1 (approximation).
#'
#' @examples
#' alpha_vals <- c(0.004, 0.005)
#' beta_vals <- c(5.0, 6.0)
#' result <- E.P.Leakage.trbb.semi.exact(
#'   alpha = alpha_vals,
#'   beta = beta_vals,
#'   cutoff = 0.02,
#'   B = 8000, b = 13,
#'   m = 5, M = 40,
#'   Delta = 0.70
#' )
#' print(result)
#'
#' @export
E.P.Leakage.trbb <- function(alpha, beta, cutoff, B, b, m, M, Delta = 1) {

  # Function to compute the incomplete Beta integral over [k, 1]
  incomplete_beta_k1 <- function(alpha, beta, k) {
    # Validate inputs
    if (any(k < 0 | k >= 1)) stop("k must be in [0, 1)")
    if (any(alpha <= 0 | beta <= 0)) stop("alpha and beta must be > 0")

    # B_{(k,1)}(alpha, beta) = B(alpha, beta) * (1 - pbeta(k, alpha, beta))
    beta_complete <- beta(alpha, beta)
    tail_prob <- 1 - pbeta(k, alpha, beta)

    return(beta_complete * tail_prob)
  }

  # Function to compute the incomplete Beta integral over [k, 1] on log scale

  incomplete_lbeta_k1 <- function(alpha, beta, k,min_log=-1e5) {
    # Validate inputs
    if (any(k < 0 | k >= 1)) stop("k must be in [0, 1)")
    if (any(alpha <= 0 | beta <= 0)) stop("alpha and beta must be > 0")

    # Compute log of B_{(k,1)}(alpha, beta)
    log_B_complete <- lbeta(alpha, beta)

    # Compute log(1 - pbeta(k, alpha, beta)) safely
    tail_prob <- 1 - pbeta(k, alpha, beta)
    log_tail_prob <- ifelse(tail_prob > 0, log(tail_prob), min_log)  # replace log(0) with floor

    #    log_tail_prob <- log1p(-pbeta(k, alpha, beta))  # log(1 - pbeta(k, ...))

    log_B_k1 <- log_B_complete + log_tail_prob

    return((log_B_k1))
  }


  bmDelta <- b * m * Delta
  E_Li <- numeric(length(alpha))
  P_Li <- numeric(length(alpha))

  for (i in seq_along(alpha)) {
    a <- alpha[i]
    b_ <- beta[i]

    # Expected Leakage: truncated beta integral
    num_trunc <- exp(incomplete_lbeta_k1(alpha = a + 1, beta = b_ + bmDelta, k = cutoff))
    denom_trunc <- exp(incomplete_lbeta_k1(alpha = a, beta = b_, k = 0))
    beta_ratio <- num_trunc / denom_trunc
    E_Li[i] <- (B - b) * M * beta_ratio

    # Probability of Leakage
    beta_mDelta <- b_ / (m * Delta)
    beta_M <- b_ / M


    if(Delta!=1.0){

      min_phi_m <- Delta * (1 - (1 - cutoff)^m)
      min_phi_M <- (1 - (1 - cutoff)^M)


      # First part: Approximation where we use beta distribution for phi_Delta
      num1 <- exp(incomplete_lbeta_k1(alpha = a, beta = beta_mDelta + b, k = min_phi_m))
      denom1 <- exp(incomplete_lbeta_k1(alpha = a, beta = beta_mDelta, k = 0))
      term1 <- num1 / denom1

      # Second part: Approximation
      # When cutoff is used, there is no difference between exact and approximate

      num2 <- exp(incomplete_lbeta_k1(alpha = a, beta = beta_M + B, k = min_phi_M))
      denom2 <- exp(incomplete_lbeta_k1(alpha = a, beta = beta_M, k = 0))
      term2 <- num2 / denom2

      # Second part: Exact
      # When cutoff is used, there is no difference between exact and approximate

      num2 <- exp(incomplete_lbeta_k1(alpha = a, beta = b_ + M*B, k = cutoff))
      denom2 <- exp(incomplete_lbeta_k1(alpha = a, beta = b_, k = 0))
      term2 <- num2 / denom2

      P_Li[i] <- term1 - term2
    } else {


      num1 <- exp(incomplete_lbeta_k1(alpha = a, beta = b_ + m*b, k = cutoff))
      denom1 <- exp(incomplete_lbeta_k1(alpha = a, beta = b_, k = 0))
      term1 <- num1 / denom1

      num2 <- exp(incomplete_lbeta_k1(alpha = a, beta = b_ + M*B, k = cutoff))
      denom2 <- exp(incomplete_lbeta_k1(alpha = a, beta = b_, k = 0))
      term2 <- num2 / denom2

      P_Li[i] <- term1 - term2
    }

  }

  return(data.frame(alpha = alpha, beta = beta, E_Li = E_Li, P_Li = P_Li))
}


#' Log-Likelihood Function for Group Testing with Imperfect Tests
#'
#' Computes the log-likelihood under a truncated Beta-Binomial model in the context of grouped sampling with imperfect sensitivity and specificity.
#' This function supports flexible parameterizations (e.g., in terms of alpha, beta, mu, rho) and includes extensions for truncated prevalence modeling,
#' numerical integration over latent prevalence values, and optional leakage estimation under imperfect sensitivity.
#' Suitable for use with optimization routines (e.g., \code{optim}) to estimate model parameters from grouped test outcome data.
#'
#' @param par A numeric vector of parameters (on log or logit scale) depending on what is fixed or estimated.
#' @param ty Observed number of positive group tests per replicate.
#' @param freq Frequency of each ty value.
#' @param theta Precision parameter. Use `Inf` for simplified form.
#' @param mu Optional; mean of the Beta prior for contamination probability.
#' @param alpha Optional; shape parameter of the Beta prior.
#' @param beta Optional; shape parameter of the Beta prior.
#' @param cutoff Truncation cutoff for the Beta distribution.
#' @param sensitivity Test sensitivity.
#' @param specificity Test specificity.
#' @param b Number of groups per replicate.
#' @param R Number of integration points for numerical approximation.
#' @param m Mean number of items per group.
#' @param M Reference sample size for leakage modeling.
#' @param deviance Logical. If TRUE, returns deviance (-2 * log-likelihood).
#' @param subtract A constant to subtract from the log-likelihood.
#' @param B Reference number of groups for leakage modeling.
#' @param mushape Logical. If TRUE, parameterize using mu and shape.
#' @param alpha.mu Logical. If TRUE, parameterize using alpha and mu.
#' @param leakage Logical. If TRUE, compute expected leakage quantities.
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
#'# For estimating log(alpha) and log(beta) for a given cutt-off greater than 0, which reflects BB model
#'mle.est.cutoff.01 <- optim(c(0, 0), loglikfn.JABES.sp.sn.TrBB,ty = ty,freq = freq,b = 13,m = 5,M = 40,theta = Inf,R = 1e4,sensitivity = 1, specificity = 1,cutoff = 0.01,deviance = FALSE,control = list(reltol = 1e-12, fnscale = -1),hessian = FALSE)
#'mle.est.cutoff.01$par
#'# For estimating log(alpha) and log(beta) for a given cutt-off of 1%, which reflects truncated BB model
#'mle.est.cutoff.01 <- optim(c(0, 0), loglikfn.JABES.sp.sn.TrBB,ty = ty,freq = freq,b = 13,m = 5,M = 40,theta = Inf,R = 1e4,sensitivity = 1, specificity = 1,cutoff = 0.01,deviance = FALSE,control = list(reltol = 1e-12, fnscale = -1),hessian = FALSE)
#'mle.est.cutoff.01$par
#'
#' @export
loglikfn.JABES.sp.sn.TrBB <- function(par,ty,freq,theta,mu,alpha,beta,cutoff,sensitivity, specificity,b,R=1000,m,M,deviance=FALSE,
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

  if(missing(sensitivity) & missing(specificity) & missing.alpha & missing(mu) & missing(beta)) {
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
      P.L <- exp(lbeta(alpha, beta / (Nbar * sensitivity) + b) - lbeta(alpha, beta / (Nbar * sensitivity))) -
        exp(lbeta(alpha, beta + B*M) - lbeta(alpha, beta))
      E.L <- (B - b) * M * exp(lbeta(alpha + 1, beta + b * Nbar * sensitivity) - lbeta(alpha, beta))
    }
  }
  if((!leakage)&deviance) return(-2*out-subtract)
  if((!leakage)&(!deviance)) return(out-subtract)
  if(leakage) return(c(deviance=-2*out,P.L=P.L,E.L=E.L))
}

