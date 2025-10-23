# --------------------------------------------------------------------------- #
# Load Required Libraries
# --------------------------------------------------------------------------- #
library(groupedBBMH)  # Custom package for grouped beta-binomial MH modeling

#-------------------------------------------------------------
# Deidentified Frozen Seafood Data
#-------------------------------------------------------------@

df.seafood <- read.csv("JASA Submission//deidentified_frozen_seafood.csv")
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

b <- 13
Nbar <- 5
ty <- seafood.data$ty
D <- length(ty)
ty2 <- x
wt <- freq
d <- length(ty2)

# Implementation of R Functions: Related to MH algorithm ------------

?MH_Sampler_BB
?summary_mcmc
?estimate_leakage
?create_mcmc_BP

G=500;R=1e3;
par <- c(log(0.005), logit(7e-04)) ; initial <- par; k <- 0.0000
MH.alpha.mu.TBB.0 <- MH_Sampler_BB(par=par,cutoff=k,
                                          initial=par,
                                          alpha=TRUE,
                                          beta=FALSE,
                                          mu=TRUE,
                                          rho=FALSE,
                                          G=G,
                                          R=R,
                                          sigma=0.2,
                                          num_chains=4,burnin=R*0.5,thin=5,
                                          seed=c(123,456,789,135),
                                          trail=b,
                                          ty2=ty2,
                                          wt=wt,
                                          group.size=Nbar,
                                          sensitivity=FALSE,
                                          specificity=FALSE,
                                          sensitivity.range=NULL,
                                          specificity.range=NULL)
summ.TBB.0 <- summary_mcmc(MH.alpha.mu.TBB.0$target.parameters$alpha_sample,
                           MH.alpha.mu.TBB.0$target.parameters$beta_sample,
                           MH.alpha.mu.TBB.0$target.parameters$mu_sample,
                          varnames = c("alpha","beta","mu"))
summ.TBB.0$summary_mcmc

G=500;R=1e3;
par <- c(log(0.005), logit(7e-04)) ; initial <- par; k <- 0.0025
MH.alpha.mu.TBB.025 <- MH_Sampler_BB(par=par,
                                     cutoff=k,
                                     initial=par,
                                     alpha=TRUE,
                                     beta=FALSE,
                                     mu=TRUE,
                                     rho=FALSE,
                                     G=G,R=R,
                                     sigma=0.2,
                                     num_chains=4,burnin=R*0.5,thin=5,
                                     seed=c(123,456,789,135),
                                     trail=b,ty2=ty2,wt=wt,group.size=Nbar,
                                     sensitivity=FALSE,
                                     specificity=FALSE,
                                     sensitivity.range=NULL,
                                     specificity.range=NULL)

# Leakage Estimate
Leakage.TBB.025 <- estimate_leakage(alpha=MH.alpha.mu.TBB.025$target.parameters$alpha_sample,
                                    beta = MH.alpha.mu.TBB.025$target.parameters$beta_sample,
                                    b=13,
                                    B=8000,
                                    m=5,
                                    M=40,sensitivity = 1.0,
                                    cutoff = 0.0025)

# Summary Statistics
summ.TBB.025 <- summary_mcmc(MH.alpha.mu.TBB.025$target.parameters$alpha_sample,
                           MH.alpha.mu.TBB.025$target.parameters$beta_sample,
                           MH.alpha.mu.TBB.025$target.parameters$mu_sample,
                           Leakage.TBB.025$E.L,
                           Leakage.TBB.025$P.L,
                           varnames = c("alpha","beta","mu","E.L","P.L"))
summ.TBB.025$summary_mcmc


# Diagnostics
mcmc_array <- create_mcmc_BP(mcmc.list.chain=MH.alpha.mu.TBB.025$parameters,
                             varnames =  c("log(alpha)", "logit(mu)"))
bayesplot::mcmc_trace(mcmc_array)
bayesplot::mcmc_dens_overlay(mcmc_array)


# Implementation of R Functions: Related to MLE method ------------

?loglik_group_bb
?loglik_group_trbb
?fit_GroupedBB
?fit_trGroupedBB

freq = c(2815,9,10,6,1,3,2,0,1,2,1,0,0,0)
ty = c(0:13)

# Fit model without leakage: Grouped BB
bb.m1 <- fit_GroupedBB(ty = ty, freq = freq, b = 13, m = 5,
                       sensitivity = 1, specificity = 1)
bb.m1 <- fit_GroupedBB(ty = ty, freq = freq, b = 13, m = 5)

# Fit model with leakage (requires B & M): Grouped BB
bb.m2 <- fit_GroupedBB(ty = ty, freq = freq, b = 13, m = 5, M = 40, B = 8000,
                       sensitivity = 1, specificity = 1, leakage = TRUE)


# Fit truncated grouped beta-binomial model with cutoff = 0 (no truncation)
trbb.m1 <- fit_trGroupedBB(ty = ty, freq = freq, b = 13, m = 5, cutoff = 0,R = 10000)
print(trbb.m1)

# Fit truncated grouped beta-binomial model with cutoff = 0.01 and imperfect sensitivity
trbb.m2 <- fit_trGroupedBB(ty, freq, b = 13, m = 5, cutoff = 0.01,
                           sensitivity = 0.8, specificity = 1,R = 10000)
print(trbb.m2)


# Estimation of Beta and Test paramaters: Application of LogLikehood function

# 5. alpha, mu, and specificity missing: : sensitivity given
MLE.BB.m5 <- optim(c(log(0.005), log(3.0), log(0.5)), loglik_group_bb,
                   ty=ty, freq=freq, b=13, theta=Inf, m=5, M=40,
                   sensitivity=1, deviance=FALSE,
                   control=list(factr=1e-12, fnscale=-1, trace=TRUE),
                   R=1e4, hessian=FALSE, method="L-BFGS-B",
                   lower=c(log(0.0001), log(0.0001), log(0.99)),
                   upper=c(Inf, Inf, 0))
MLE.BB.m5$mu <- exp(MLE.BB.m5$par[1]) / sum(exp(MLE.BB.m5$par))
MLE.BB.m5$alpha <- exp(MLE.BB.m5$par[1])
MLE.BB.m5$beta <- exp(MLE.BB.m5$par[2])
MLE.BB.m5$specificity <- exp(MLE.BB.m5$par[3])


# 6. alpha, mu, sensitivity, and specificity are all missing
MLE.BB.m6 <- optim(c(0,0,0,0), loglik_group_bb,
                      ty=ty, freq=freq, b=13, theta=Inf, m=5, M=40,
                      deviance=FALSE,
                      control=list(factr=1e-12, fnscale=-1, trace=TRUE),
                      R=1e4, hessian=FALSE, method="L-BFGS-B",
                      lower=c(log(0.0001), log(0.0001), log(0.5), log(0.99)),
                      upper=c(Inf, Inf, 0, 0))
MLE.BB.m6$mu <- exp(MLE.BB.m6$par[1]) / sum(exp(MLE.BB.m6$par))
MLE.BB.m6$alpha <- exp(MLE.BB.m6$par[1])
MLE.BB.m6$beta <- exp(MLE.BB.m6$par[2])
MLE.BB.m6$sensitivity <- exp(MLE.BB.m6$par[3])
MLE.BB.m6$specificity <- exp(MLE.BB.m6$par[4])
# 6. alpha, mu, sensitivity, and specificity are all missing
MLE.BB.m6 <- optim(c(0,0,0,0), loglik_group_bb,
                   ty=ty, freq=freq, b=13, theta=Inf, m=5, M=40,
                   deviance=FALSE,
                   control=list(factr=1e-12, fnscale=-1, trace=TRUE),
                   R=1e4, hessian=FALSE, method="L-BFGS-B",
                   lower=c(log(0.0001), log(0.0001), log(0.5), log(0.99)),
                   upper=c(Inf, Inf, 0, 0))
MLE.BB.m6$mu <- exp(MLE.BB.m6$par[1]) / sum(exp(MLE.BB.m6$par))
MLE.BB.m6$alpha <- exp(MLE.BB.m6$par[1])
MLE.BB.m6$beta <- exp(MLE.BB.m6$par[2])
MLE.BB.m6$sensitivity <- exp(MLE.BB.m6$par[3])
MLE.BB.m6$specificity <- exp(MLE.BB.m6$par[4])



# For estimating log(beta) for a given alpha and cutoff greater than 0, which reflects truncated BB model
beta.est.cutoff.01 <- optim(c(0),loglik_group_trbb,ty=ty,freq=freq,b=13,m=5,theta=Inf,R=1e4,alpha=0.005,cutoff = 0.01,sensitivity=1,specificity=1,
                            deviance=FALSE,method="Brent", control=list(reltol=1e-12,fnscale=-1), hessian=FALSE,lower = c(log(0.05)),upper = c(log(100)))
exp(beta.est.cutoff.01$par)

