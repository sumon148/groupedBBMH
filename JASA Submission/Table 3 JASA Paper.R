
# --------------------------------------------------------------------------- #
# Load Required Libraries
# --------------------------------------------------------------------------- #
library(groupedBBMH)  # Custom package for grouped beta-binomial MH modeling
library(dplyr)        # Data manipulation
library(xtable)       # LaTeX table generation
library(MASS)         # For multivariate normal distribution functions
library(extraDistr)   # Additional distributions (e.g., beta-binomial, truncated distributions)
library(reshape2)     # Data reshaping (e.g., melt and dcast for wide/long format)
library(bayesplot)    # Diagnostic and visualization tools for MCMC
library(coda)         # MCMC diagnostics and summaries
library(brms)         # Bayesian regression models using Stan backend
library(loo)
library(ggplot2)      # General plotting
library(ggpubr)       # Publication-ready plots and arrangement tools

#-------------------------------------------------------------@
# Table 3: JASA Paper
# Using seafood data
# Exact BB Model: Frequentist method using MLE and PLLF
# Exact BB Model: Bayesian method MH Algorithm
# ML and MH: Both are based on Generalized BB model with cutoff, k
# Approximate BB model: Bayesian method using BRMS package
# Target Parameters: alpha, beta, mu, E.L, P.l
# keep Format of Old Tables First
# Then Modify them for paper
# Relevant Supplementary Table will be from here
#-------------------------------------------------------------@

#-------------------------------------------------------------
# Deidentified Frozen Seafood Data
#-------------------------------------------------------------

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

#-----------------------------------------@
# define data for MH algorithm
b <- 13
Nbar <- 5
ty <- seafood.data$ty
D <- length(ty)
ty2 <- x
wt <- freq
d <- length(ty2)
#-----------------------------------------@


# Dispersion Index
dispersion.index <- function(y) {
  DI <- var(y,na.rm=T) / mean(y,na.rm=T)
  return(DI)
}

var(seafood.data$ty) / mean(seafood.data$ty)
mean(seafood.data$ty[seafood.data$ty>0])
median(seafood.data$ty[seafood.data$ty>0])
var(seafood.data$ty[seafood.data$ty>0]) / mean(seafood.data$ty[seafood.data$ty>0])
dispersion.index(seafood.data$ty[seafood.data$ty>0])

# Exact BB Model: MLE --------------------------------

# ----------------------------------------------------------------------#
# MLE Estimate: alpha, beta, mu, E.L and P.L - Perfect test accuracy case
# ----------------------------------------------------------------------#

o.theta.inf.seafood <- optim(c(0,0),loglik_group_trbb ,ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,
                           theta=Inf,m=5,M=40,deviance=FALSE,control=list(reltol=1e-12,fnscale=-1),
                           R=1e4,hessian=FALSE,sensitivity=1,specificity=1,cutoff = 0)
o.theta.inf.seafood$mu <- exp(o.theta.inf.seafood$par[1])/sum(exp(o.theta.inf.seafood$par))
o.theta.inf.seafood$alpha <- exp(o.theta.inf.seafood$par[1])
o.theta.inf.seafood$beta <- exp(o.theta.inf.seafood$par[2])


temp <- loglik_group_trbb(o.theta.inf.seafood$par,ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,B=8000,
                             theta=Inf,m=5,M=40,deviance=FALSE,leakage=TRUE,R=1e4,
                             sensitivity=1,specificity=1,cutoff = 0)
o.theta.inf.seafood$P.L <- temp["P.L"]
o.theta.inf.seafood$E.L <- temp["E.L"]

# Check leakage calculation
E.P.Leakge.sn.1 <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                    alpha=o.theta.inf.seafood$alpha,beta=o.theta.inf.seafood$beta,cutoff = 0,
                                    sensitivity = 1)
E.P.Leakge.sn.1$mu <- E.P.Leakge.sn.1$alpha/(E.P.Leakge.sn.1$alpha+E.P.Leakge.sn.1$beta)

MLE.perfect.case <- cbind(E.P.Leakge.sn.1$alpha,
                          E.P.Leakge.sn.1$beta,
                          E.P.Leakge.sn.1$mu,
                          E.P.Leakge.sn.1$E.L,
                          E.P.Leakge.sn.1$P.L)
colnames(MLE.perfect.case) <- c("alpha","beta","mu","E.L","P.L")
rownames(MLE.perfect.case) <- c("Perfect Sensitivity")
MLE.perfect.case
#                         alpha     beta           mu      E.L       P.L
# Perfect Sensitivity 0.005125782 6.883926 0.0007440476 22.49933 0.04166347



# ----------------------------------------------------------------------#
# MLE Estimate: alpha, beta, mu, E.L and P.L - Imperfect test  (Delta=0.80) accuracy case
# ----------------------------------------------------------------------#

# alpha and mu missing
o.theta.inf.seafood.sn.80 <- optim(c(0,0),loglik_group_trbb,ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,
                                 theta=Inf,m=5,M=40,deviance=FALSE,control=list(reltol=1e-12,fnscale=-1),
                                 R=1e4,hessian=FALSE,sensitivity=0.80,specificity=1,cutoff = 0)
o.theta.inf.seafood.sn.80$mu <- exp(o.theta.inf.seafood.sn.80$par[1])/sum(exp(o.theta.inf.seafood.sn.80$par))
o.theta.inf.seafood.sn.80$alpha <- exp(o.theta.inf.seafood.sn.80$par[1])
o.theta.inf.seafood.sn.80$beta <- exp(o.theta.inf.seafood.sn.80$par[2])

temp <- loglik_group_trbb(o.theta.inf.seafood.sn.80$par,ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,B=8000,
                             theta=Inf,m=5,M=40,deviance=FALSE,leakage=TRUE,R=1e4,sensitivity=0.80,specificity=1,cutoff = 0)
o.theta.inf.seafood.sn.80$P.L <- temp["P.L"]
o.theta.inf.seafood.sn.80$E.L <- temp["E.L"]

E.P.Leakge.sn.0.8 <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                    alpha=o.theta.inf.seafood.sn.80$alpha,beta=o.theta.inf.seafood.sn.80$beta,cutoff = 0,
                                    sensitivity = 0.8)
E.P.Leakge.sn.0.8$mu <- E.P.Leakge.sn.0.8$alpha/(E.P.Leakge.sn.0.8$alpha+E.P.Leakge.sn.0.8$beta)

MLE.imperfect.case.sn.80 <- cbind(E.P.Leakge.sn.0.8$alpha,
                                  E.P.Leakge.sn.0.8$beta,
                                  E.P.Leakge.sn.0.8$mu,
                                  E.P.Leakge.sn.0.8$E.L,
                                  E.P.Leakge.sn.0.8$P.L)
colnames(MLE.imperfect.case.sn.80) <- c("alpha","beta","mu","E.L","P.L")
rownames(MLE.imperfect.case.sn.80) <- c("Imperfect Sensitivity")
# MLE.imperfect.case.sn.80
#                         alpha     beta          mu      E.L       P.L
# Imperfect Sensitivity 0.004679024 4.398965 0.001062535 26.17504 0.03747001

# ----------------------------------------------------------------------
# Profile Likelihood Confidence interval - alpha : Perfect and Imperfect test accuracy
# Profile likelihood: Beta and Leakage Parameters -------------#
# alpha and mu for varying sensitivity (0.70,0.75,0.80, 0.90,0.95,0.99,1)
# Perfect specificity
# Results are already saved. Don't run the code
# ----------------------------------------------------------------------

# sensitivity  <- c(0.70,0.75,0.80,0.90,0.95,0.99,1)
sensitivity  <- c(0.80)
specificity  <- 1
alpha  <- seq(0.002,0.008,by=0.00001) # 0.00001
df.pllf.sn.sp.alpha <- expand.grid(sensitivity,specificity,alpha)
dim(df.pllf.sn.sp.alpha)

colnames(df.pllf.sn.sp.alpha) <- c("sensitivity","specificity","alpha")
df.pllf.sn.sp.alpha$logL <- rep(NA,dim(df.pllf.sn.sp.alpha)[1])
# df.pllf.sn.sp.alpha$alpha <- rep(NA,dim(df.pllf.sn.sp.alpha)[1])
df.pllf.sn.sp.alpha$beta <- rep(NA,dim(df.pllf.sn.sp.alpha)[1])
df.pllf.sn.sp.alpha$mu <- rep(NA,dim(df.pllf.sn.sp.alpha)[1])
df.pllf.sn.sp.alpha$converge <- rep(NA,dim(df.pllf.sn.sp.alpha)[1])
df.pllf.sn.sp.alpha$P.L <- rep(NA,dim(df.pllf.sn.sp.alpha)[1])
df.pllf.sn.sp.alpha$E.L <- rep(NA,dim(df.pllf.sn.sp.alpha)[1])


for (i in 1:dim(df.pllf.sn.sp.alpha)[1]){

  skip_to_next <- FALSE

  tryCatch(

    {o.theta.inf.seafood.sp <- optim(c(0),loglik_group_trbb,ty=c(summ.seafood.data$ty),freq=c(summ.seafood.data$freq),b=13,m=5,
                                   theta=Inf,R=1e4,alpha=df.pllf.sn.sp.alpha$alpha[i],
                                   sensitivity=df.pllf.sn.sp.alpha$sensitivity[i],
                                   specificity=df.pllf.sn.sp.alpha$specificity[i],
                                   cutoff=0,
                                   deviance=FALSE,
                                   method="Brent", control=list(reltol=1e-12,fnscale=-1),
                                   hessian=FALSE,
                                   lower = c(log(0.05)),upper = c(log(100)))

    temp <- loglik_group_trbb(c(log(df.pllf.sn.sp.alpha$alpha[i]),o.theta.inf.seafood.sp$par),ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,B=8000,
                                 theta=Inf,m=5,M=40,leakage=TRUE,R=1e4,
                              sensitivity=df.pllf.sn.sp.alpha$sensitivity[i],
                              specificity=df.pllf.sn.sp.alpha$specificity[i],
                              cutoff=0)
    o.theta.inf.seafood.sp$P.L <- temp["P.L"]
    o.theta.inf.seafood.sp$E.L <- temp["E.L"]



    df.pllf.sn.sp.alpha$beta[i] <- exp(o.theta.inf.seafood.sp$par[1])
    df.pllf.sn.sp.alpha$mu[i] <- df.pllf.sn.sp.alpha$alpha[i]/sum(df.pllf.sn.sp.alpha$alpha[i], df.pllf.sn.sp.alpha$beta[i])

    df.pllf.sn.sp.alpha$logL[i] <- o.theta.inf.seafood.sp$value
    df.pllf.sn.sp.alpha$converge[i] <- o.theta.inf.seafood.sp$convergence
    df.pllf.sn.sp.alpha$P.L[i] <- o.theta.inf.seafood.sp$P.L
    df.pllf.sn.sp.alpha$E.L[i] <- o.theta.inf.seafood.sp$E.L
    },
    error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }


}


df.pllf.sn.sp.alpha.prawn <- df.pllf.sn.sp.alpha

save(df.pllf.sn.sp.alpha.prawn,file="df.pllf.sn.sp.alpha.prawn.updated.rdata")

# save(df.pllf.sn.sp.alpha.prawn,file="df.pllf.sn.sp.alpha.prawn.rdata")
# save(df.pllf.sn.sp.alpha.prawn,file="df.pllf.sn.sp.alpha.prawn.80.rdata")


# Leakage parameter estimates: Sensitivity 1.0 and 0.8

load("df.pllf.sn.sp.alpha.prawn.updated.rdata")

# ----------------------------------------------------------------------
# Profile Likelihood Confidence interval - mu : Perfect and Imperfect test accuracy
# Profile likelihood: Beta and Leakage Parameters
# alpha and mu for varying sensitivity (0.70,0.75,0.80, 0.90,0.95,0.99,1)
# Under perfect specificity --------------#
# Results are already saved. Don't run the code
# ----------------------------------------------------------------------

sensitivity  <- c(0.70,0.75,0.80, 0.90,0.95,0.99,1)
# sensitivity  <- c(0.80)
specificity  <- 1
mu  <- seq(0.0003,0.004,by=0.00001) # 00001
df.pllf.sn.sp.mu <- expand.grid(sensitivity,specificity,mu)
dim(df.pllf.sn.sp.mu)

colnames(df.pllf.sn.sp.mu) <- c("sensitivity","specificity","mu")
df.pllf.sn.sp.mu$logL <- rep(NA,dim(df.pllf.sn.sp.mu)[1])
df.pllf.sn.sp.mu$alpha <- rep(NA,dim(df.pllf.sn.sp.mu)[1])
df.pllf.sn.sp.mu$beta <- rep(NA,dim(df.pllf.sn.sp.mu)[1])
# df.pllf.sn.sp.mu$mu <- rep(NA,dim(df.pllf.sn.sp.mu)[1])
df.pllf.sn.sp.mu$converge <- rep(NA,dim(df.pllf.sn.sp.mu)[1])
df.pllf.sn.sp.mu$P.L <- rep(NA,dim(df.pllf.sn.sp.mu)[1])
df.pllf.sn.sp.mu$E.L <- rep(NA,dim(df.pllf.sn.sp.mu)[1])



for (i in 1:dim(df.pllf.sn.sp.mu)[1]){

  skip_to_next <- FALSE

  tryCatch(

    {o.theta.inf.seafood.sp <- optim(c(0),loglik_group_trbb,ty=c(summ.seafood.data$ty),freq=c(summ.seafood.data$freq),b=13,m=5,
                                   theta=Inf,R=1e4,mu=df.pllf.sn.sp.mu$mu[i], sensitivity=df.pllf.sn.sp.mu$sensitivity[i],specificity=df.pllf.sn.sp.mu$specificity[i],
                                   deviance=FALSE,method="Brent", control=list(reltol=1e-12,fnscale=-1), hessian=FALSE,lower = c(log(0.0005)),upper = c(log(100)))

    o.theta.inf.seafood.sp$alpha <- exp(o.theta.inf.seafood.sp$par)
    o.theta.inf.seafood.sp$beta <- (1-df.pllf.sn.sp.mu$mu[i])/df.pllf.sn.sp.mu$mu[i]*o.theta.inf.seafood.sp$alpha
    temp <- loglik_group_trbb(c(log(o.theta.inf.seafood.sp$alpha),log(o.theta.inf.seafood.sp$beta)),ty=summ.seafood.data$ty,freq=summ.seafood.data$freq,b=13,B=8000,
                                 theta=Inf,m=5,M=40,leakage=TRUE,R=1e4,sensitivity=df.pllf.sn.sp.mu$sensitivity[i],specificity=df.pllf.sn.sp.mu$specificity[i])
    o.theta.inf.seafood.sp$P.L <- temp["P.L"]
    o.theta.inf.seafood.sp$E.L <- temp["E.L"]

    df.pllf.sn.sp.mu$logL[i] <- o.theta.inf.seafood.sp$value
    df.pllf.sn.sp.mu$alpha[i] <- o.theta.inf.seafood.sp$alpha
    df.pllf.sn.sp.mu$beta[i] <- o.theta.inf.seafood.sp$beta
    df.pllf.sn.sp.mu$converge[i] <- o.theta.inf.seafood.sp$convergence
    df.pllf.sn.sp.mu$P.L[i] <- o.theta.inf.seafood.sp$P.L
    df.pllf.sn.sp.mu$E.L[i] <- o.theta.inf.seafood.sp$E.L},

    error = function(e) { skip_to_next <<- TRUE})

  if(skip_to_next) { next }


}


df.pllf.sn.sp.mu.prawn <- df.pllf.sn.sp.mu

# save(df.pllf.sn.sp.mu.prawn,file="df.pllf.sn.sp.mu.prawn.updated.rdata")
# save(df.pllf.sn.sp.mu.prawn,file="df.pllf.sn.sp.mu.prawn.rdata")
# save(df.pllf.sn.sp.mu.prawn,file="df.pllf.sn.sp.mu.prawn.80.rdata")
save(df.pllf.sn.sp.mu.prawn,file="df.pllf.sn.sp.mu.prawn.updated.rdata")


#-------------------------------------------
# Load PLLF for alpha and mu
# PLLF for alpha and mu provide PLLF of beta
# PLLF for alpha provide PLLF for E.L and P.L
#-------------------------------------------


load("df.pllf.sn.sp.alpha.prawn.updated.rdata")
load("df.pllf.sn.sp.mu.prawn.updated.rdata")


df.pllf.sn.sp.alpha.prawn$P.L[df.pllf.sn.sp.alpha.prawn$sensitivity==1.0] <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                    alpha=df.pllf.sn.sp.alpha.prawn$alpha[df.pllf.sn.sp.alpha.prawn$sensitivity==1.0],beta=df.pllf.sn.sp.alpha.prawn$beta[df.pllf.sn.sp.alpha.prawn$sensitivity==1.0],cutoff = 0,
                                    Delta=1,approximate=FALSE)$P_Li
df.pllf.sn.sp.alpha.prawn$E.L[df.pllf.sn.sp.alpha.prawn$sensitivity==1.0] <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                   alpha=df.pllf.sn.sp.alpha.prawn$alpha[df.pllf.sn.sp.alpha.prawn$sensitivity==1.0],beta=df.pllf.sn.sp.alpha.prawn$beta[df.pllf.sn.sp.alpha.prawn$sensitivity==1.0],cutoff = 0,
                                   Delta=1,approximate=FALSE)$E_Li
df.pllf.sn.sp.alpha.prawn$P.L[df.pllf.sn.sp.alpha.prawn$sensitivity==0.8] <-  estimate_leakage(B=8000,b=13,M=40,m=5,
                                   alpha=df.pllf.sn.sp.alpha.prawn$alpha[df.pllf.sn.sp.alpha.prawn$sensitivity==0.8],
                                   beta=df.pllf.sn.sp.alpha.prawn$beta[df.pllf.sn.sp.alpha.prawn$sensitivity==0.8],
                                   cutoff = 0,
                                   Delta=0.8,approximate=TRUE)$P_Li
df.pllf.sn.sp.alpha.prawn$E.L[df.pllf.sn.sp.alpha.prawn$sensitivity==0.8] <-  estimate_leakage(B=8000,b=13,M=40,m=5,
                                   alpha=df.pllf.sn.sp.alpha.prawn$alpha[df.pllf.sn.sp.alpha.prawn$sensitivity==0.8],
                                   beta=df.pllf.sn.sp.alpha.prawn$beta[df.pllf.sn.sp.alpha.prawn$sensitivity==0.8],
                                   cutoff = 0,
                                   Delta=0.8,approximate=TRUE)$E_Li


pllf.CI.alpha.sn.1 <- pllfCIestimate(df.pllf.sn.sp.alpha.prawn$alpha[df.pllf.sn.sp.alpha.prawn$sensitivity==1],
                                       logL=df.pllf.sn.sp.alpha.prawn$logL[df.pllf.sn.sp.alpha.prawn$sensitivity==1],conf_level = 0.95)
pllf.CI.mu.sn.1 <- pllfCIestimate(df.pllf.sn.sp.mu.prawn$mu[df.pllf.sn.sp.mu.prawn$sensitivity==1],
                                    logL=df.pllf.sn.sp.mu.prawn$logL[df.pllf.sn.sp.mu.prawn$sensitivity==1])
pllf.CI.beta.sn.1 <- pllfCIestimate(df.pllf.sn.sp.mu.prawn$beta[df.pllf.sn.sp.mu.prawn$sensitivity==1],
                                      logL=df.pllf.sn.sp.mu.prawn$logL[df.pllf.sn.sp.mu.prawn$sensitivity==1],conf_level = 0.95)
pllf.CI.E.L.sn.1 <- pllfCIestimate(df.pllf.sn.sp.alpha.prawn$E.L[df.pllf.sn.sp.alpha.prawn$sensitivity==1],
                                     logL=df.pllf.sn.sp.alpha.prawn$logL[df.pllf.sn.sp.alpha.prawn$sensitivity==1],conf_level = 0.95)
pllf.CI.P.L.sn.1 <- pllfCIestimate(df.pllf.sn.sp.alpha.prawn$P.L[df.pllf.sn.sp.alpha.prawn$sensitivity==1],
                                     logL=df.pllf.sn.sp.alpha.prawn$logL[df.pllf.sn.sp.alpha.prawn$sensitivity==1],conf_level = 0.95)
pllf.CI.alpha.sn.0.80 <- pllfCIestimate(df.pllf.sn.sp.alpha.prawn$alpha[df.pllf.sn.sp.alpha.prawn$sensitivity==0.80],
                                          logL=df.pllf.sn.sp.alpha.prawn$logL[df.pllf.sn.sp.alpha.prawn$sensitivity==0.80],conf_level = 0.95)
pllf.CI.mu.sn.0.80 <- pllfCIestimate(df.pllf.sn.sp.mu.prawn$mu[df.pllf.sn.sp.mu.prawn$sensitivity==0.80],
                                       logL=df.pllf.sn.sp.mu.prawn$logL[df.pllf.sn.sp.mu.prawn$sensitivity==0.80])
pllf.CI.beta.sn.0.80 <- pllfCIestimate(df.pllf.sn.sp.mu.prawn$beta[df.pllf.sn.sp.mu.prawn$sensitivity==0.80],
                                         logL=df.pllf.sn.sp.mu.prawn$logL[df.pllf.sn.sp.mu.prawn$sensitivity==0.80],conf_level = 0.95)
pllf.CI.E.L.sn.0.80 <- pllfCIestimate(df.pllf.sn.sp.alpha.prawn$E.L[df.pllf.sn.sp.alpha.prawn$sensitivity==0.80],
                                        logL=df.pllf.sn.sp.alpha.prawn$logL[df.pllf.sn.sp.alpha.prawn$sensitivity==0.80],conf_level = 0.95)
pllf.CI.P.L.sn.0.80 <- pllfCIestimate(df.pllf.sn.sp.alpha.prawn$P.L[df.pllf.sn.sp.alpha.prawn$sensitivity==0.80],
                                        logL=df.pllf.sn.sp.alpha.prawn$logL[df.pllf.sn.sp.alpha.prawn$sensitivity==0.80],conf_level = 0.95)


#=====================================================#
# Results from MLE and PLLF: Perfect and imperfect Sensitivity
#=====================================================#

alpha_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"alpha"],pllf.CI.alpha.sn.1$mle.ll,pllf.CI.alpha.sn.1$mle.ul)
alpha_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"alpha"],pllf.CI.alpha.sn.0.80$mle.ll,pllf.CI.alpha.sn.0.80$mle.ul)

beta_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"beta"],pllf.CI.beta.sn.1$mle.ll,pllf.CI.beta.sn.1$mle.ul)
beta_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"beta"],pllf.CI.beta.sn.0.80$mle.ll,pllf.CI.beta.sn.0.80$mle.ul)

mu_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"mu"],pllf.CI.mu.sn.1$mle.ll,pllf.CI.mu.sn.1$mle.ul)
mu_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"mu"],pllf.CI.mu.sn.0.80$mle.ll,pllf.CI.mu.sn.0.80$mle.ul)

E.L_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"E.L"],pllf.CI.E.L.sn.1$mle.ll,pllf.CI.E.L.sn.1$mle.ul)
E.L_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"E.L"],pllf.CI.E.L.sn.0.80$mle.ll,pllf.CI.E.L.sn.0.80$mle.ul)

P.L_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"P.L"],pllf.CI.P.L.sn.1$mle.ll,pllf.CI.P.L.sn.1$mle.ul)
P.L_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"P.L"],pllf.CI.P.L.sn.0.80$mle.ll,pllf.CI.P.L.sn.0.80$mle.ul)






# ---------------------------------------------------------------
# Approx. BB Model: BRMS assuming Perfect and imperfect Sensitivity
# ---------------------------------------------------------------


m1.prawn.a <- brm(ty | trials(13) ~ 1,
                  family = beta_binomial(link = "logit", link_phi = "log"), # binomial("identity") would be more straightforward
                  data = seafood.data,
                  chains = 4, # nb of chains
                  iter = 6000, # nb of iterations, including burnin
                  warmup = 1000, # burnin
                  thin = 5) # thinning
m1.prawn.a$prior

save(m1.prawn.a,file="m1.prawn.a.Rdata")

# ---------------------------------------------------------------
# Results from Approximate model : HMC Estimator
# ---------------------------------------------------------------

load("m1.prawn.a.Rdata")

summary(m1.prawn.a,waic=TRUE)
loo.m1.prawn.a <- loo(m1.prawn.a)
loo.m1.prawn.a$estimates

## predicted responses
pp <- predict(m1.prawn.a,summary = TRUE,probs = c(0.025, 0.975))

# extract posterior draws in an array format
(draws_fit_prawn_a <- brms::as_draws_array(m1.prawn.a))
# extract two main parameters: "b_Intercept" and "phi" (theta in our notation)
b_Intercept_draws_a <- brms::as_draws_array(m1.prawn.a, variable = "b_Intercept")
phi_draws_prawn_a <- brms::as_draws_array(m1.prawn.a, variable = "phi")

# Generate P values directly from fitted model
p_draws_alt <- brms::posterior_linpred(m1.prawn.a,  transform = TRUE) # Probabilities for all batches
# Batch specific probability
p_batch <- apply(p_draws_alt,2,mean)
p_batch_median <- apply(p_draws_alt,2,median)

# We can use transformation function as well
p_draws_prawn_a <- posterior::mutate_variables(b_Intercept_draws_a, p = (1-1/(1+exp(b_Intercept))))
p_draws_prawn_a <- p_draws_prawn_a[,,2]
p_draws_vector <- as.vector(p_draws_prawn_a)
summ_p <- summary(p_draws_vector)
CI_p <- quantile(p_draws_vector,probs = c(0.0275,0.975))

summary(p_draws_prawn_a)
summary(b_Intercept_draws_a)

# Theta parameters
summ_theta_prawn <- summary(as.vector(phi_draws_prawn_a))
CI_theta_prawn <- quantile(as.vector(phi_draws_prawn_a),probs = c(0.0275,0.975))

# Shape 1 and Shape 2
shape1_draws_prawn <- p_draws_prawn_a*phi_draws_prawn_a
summ_alpha_prawn <- summary(shape1_draws_prawn)
CI_alpha_prawn <- quantile(as.vector(shape1_draws_prawn),probs = c(0.0275,0.975))
alpha_estimate_prawn <- c(as.numeric(summ_alpha_prawn[,"mean"]),as.vector(CI_alpha_prawn))

shape2_draws_prawn <- (1-p_draws_prawn_a)*phi_draws_prawn_a
summ_beta_prawn <- summary(shape2_draws_prawn)
CI_beta_prawn <- quantile(as.vector(shape2_draws_prawn),probs = c(0.0275,0.975))
beta_estimate_prawn <- c(as.numeric(summ_beta_prawn[,"mean"]),as.vector(CI_beta_prawn))


# Beta: Perfect and Imperfect Sensitivity (Delta=0.70)  ---------- #

shape2_draws_prawn_Delta_1.0 <- shape2_draws_prawn*5*1.0
summ_beta_prawn_Delta_1.0 <- summary(shape2_draws_prawn_Delta_1.0)
CI_beta_prawn_Delta_1.0 <- quantile(as.vector(shape2_draws_prawn_Delta_1.0),probs = c(0.0275,0.975))
beta_estimate_prawn_Delta_1.0 <- c(as.numeric(summ_beta_prawn_Delta_1.0[,"mean"]),as.vector(CI_beta_prawn_Delta_1.0))

shape2_draws_prawn_Delta_0.80 <- shape2_draws_prawn*5*0.80
summ_beta_prawn_Delta_0.80 <- summary(shape2_draws_prawn_Delta_0.80)
CI_beta_prawn_Delta_0.80 <- quantile(as.vector(shape2_draws_prawn_Delta_0.80),probs = c(0.0275,0.975))
beta_estimate_prawn_Delta_0.80 <- c(as.numeric(summ_beta_prawn_Delta_0.80[,"mean"]),as.vector(CI_beta_prawn_Delta_0.80))

mu_draws_prawn_Delta_0.80 <- shape1_draws_prawn/(shape1_draws_prawn+shape2_draws_prawn_Delta_0.80)
summ_mu_prawn_Delta_0.80 <- summary(mu_draws_prawn_Delta_0.80)
CI_mu_prawn_Delta_0.80 <- quantile(as.vector(mu_draws_prawn_Delta_0.80),probs = c(0.0275,0.975))
mu_estimate_prawn_Delta_0.80 <- c(as.numeric(summ_mu_prawn_Delta_0.80[,"mean"]),as.vector(CI_mu_prawn_Delta_0.80))

mu_draws_prawn_Delta_1.0 <- shape1_draws_prawn/(shape1_draws_prawn+shape2_draws_prawn_Delta_1.0)
summ_mu_prawn_Delta_1.0 <- summary(mu_draws_prawn_Delta_1.0)
CI_mu_prawn_Delta_1.0 <- quantile(as.vector(mu_draws_prawn_Delta_1.0),probs = c(0.0275,0.975))
mu_estimate_prawn_Delta_1.0 <- c(as.numeric(summ_mu_prawn_Delta_1.0[,"mean"]),as.vector(CI_mu_prawn_Delta_1.0))

# Calculation of E.L and P.L ----------#
E.L.P.L.Draws.Lambda.0.80.prawn <- estimate_leakage(B=8000,b=13,M=40,m=5,alpha=shape1_draws_prawn,beta=shape2_draws_prawn_Delta_0.80,sensitivity=0.80)
summ.E.L.Draws.Lambda.0.80.prawn <- summary(E.L.P.L.Draws.Lambda.0.80.prawn$E.L)
summ.P.L.Draws.Lambda.0.80.prawn <- summary(E.L.P.L.Draws.Lambda.0.80.prawn$P.L)
CI.E.L.Draws.Lambda.0.80.prawn <- quantile(as.vector(E.L.P.L.Draws.Lambda.0.80.prawn$E.L),probs = c(0.0275,0.975))
CI.P.L.Draws.Lambda.0.80.prawn <- quantile(as.vector(E.L.P.L.Draws.Lambda.0.80.prawn$P.L),probs = c(0.0275,0.975))
E.L_estimate_prawn_Delta_0.80 <- c(as.numeric(summ.E.L.Draws.Lambda.0.80.prawn[,"mean"]),as.vector(CI.E.L.Draws.Lambda.0.80.prawn))
P.L_estimate_prawn_Delta_0.80 <- c(as.numeric(summ.P.L.Draws.Lambda.0.80.prawn[,"mean"]),as.vector(CI.P.L.Draws.Lambda.0.80.prawn))


E.L.P.L.Draws.Lambda.1.prawn <- estimate_leakage(B=8000,b=13,M=40,m=5,alpha=shape1_draws_prawn,beta=shape2_draws_prawn_Delta_1.0,sensitivity=1)  # If no clustering but 5 seeds are pooled to test
summ.E.L.Draws.Lambda.1.prawn <- summary(E.L.P.L.Draws.Lambda.1.prawn$E.L)
summ.P.L.Draws.Lambda.1.prawn <- summary(E.L.P.L.Draws.Lambda.1.prawn$P.L)
CI.E.L.Draws.Lambda.1.prawn <- quantile(as.vector(E.L.P.L.Draws.Lambda.1.prawn$E.L),probs = c(0.0275,0.975))
CI.P.L.Draws.Lambda.1.prawn <- quantile(as.vector(E.L.P.L.Draws.Lambda.1.prawn$P.L),probs = c(0.0275,0.975))
E.L_estimate_prawn_Delta_1 <- c(as.numeric(summ.E.L.Draws.Lambda.1.prawn[,"mean"]),as.vector(CI.E.L.Draws.Lambda.1.prawn))
P.L_estimate_prawn_Delta_1 <- c(as.numeric(summ.P.L.Draws.Lambda.1.prawn[,"mean"]),as.vector(CI.P.L.Draws.Lambda.1.prawn))

#=====================================================#
# Results from BRMS: Perfect and imperfect Sensitivity
#=====================================================#

alpha_estimate_BRMS <- alpha_estimate_prawn
beta_start_estimate_BRMS <- beta_estimate_prawn   # beta_star=beta/(m*sensitivity)

beta_estimate_SN.1_BRMS <-beta_estimate_prawn_Delta_1.0 # beta=beta_star*m*sensitivity
beta_estimate_SN.0.80_BRMS <-beta_estimate_prawn_Delta_0.80 # beta=beta_star*m*sensitivity

mu_estimate_SN.1_BRMS <-mu_estimate_prawn_Delta_1.0
mu_estimate_SN.0.80_BRMS <-mu_estimate_prawn_Delta_0.80

E.L_estimate_SN.1_BRMS <-E.L_estimate_prawn_Delta_1
E.L_estimate_SN.0.80_BRMS <-E.L_estimate_prawn_Delta_0.80

P.L_estimate_SN.1_BRMS <-P.L_estimate_prawn_Delta_1
P.L_estimate_SN.0.80_BRMS <-P.L_estimate_prawn_Delta_0.80




# Exact BB Model: MH algorithm assuming Perfect and imperfect Sensitivity --------------------------------

#-----------------------------------------@
# define data
b <- 13
Nbar <- 5
ty <- seafood.data$ty
D <- length(ty)
ty2 <- x
wt <- freq
d <- length(ty2)
#------------------------------------------@

#-------------------------------------------------------------------------#
# Model under Perfect Sensitivity and Perfect specificity
# Sensitivity: 1
# Specificity: 1

G=3000;R=50e3;
par <- c(log(0.005), logit(7e-04)) ; initial <- par;
MH.alpha.mu.sigma.20.Perfect.Prawn <- MH_Sampler_BB(par=par,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                               sigma=0.2,num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=FALSE,specificity=FALSE,
                                               sensitivity.range=NULL,specificity.range=NULL)
save(MH.alpha.mu.sigma.20.Perfect.Prawn,file="MH.alpha.mu.sigma.20.Perfect.Prawn.Rdata")

#-------------------------------------------------------------------------#
# Model under known imperfect Sensitivity but perfect specificity
# Sensitivity: 0.80
# Specificity: 1

G=3000;R=50e3;
par <- c(log(0.005), logit(7e-04), logit(0.75)) ; initial <- par;
MH.alpha.mu.sigma.20.known.sn.80.Prawn <- MH_Sampler_BB(par=par,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                           sigma=c(0.2,0.2),num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=TRUE,specificity=FALSE,
                                                           sensitivity.range=c(0.80,0.80),specificity.range=NULL)
save(MH.alpha.mu.sigma.20.known.sn.80.Prawn,file="MH.alpha.mu.sigma.20.known.sn.80.Prawn.Rdata")


#-------------------------------------------------------------------------#
# Model under (unknown) imperfect Sensitivity  but perfect specificity
# Sensitivity: U(0.75, 0.85)
# Specificity: 1

G=3000;R=50e3;
par <- c(log(0.005), logit(7e-04), logit(0.75)) ; initial <- par;
MH.alpha.mu.sigma.20.unknown.sn.80.Prawn <- MH_Sampler_BB(par=par,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                          sigma=c(0.2,0.2,2.0),num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=TRUE,specificity=FALSE,
                                                          sensitivity.range=c(0.75,0.85),specificity.range=NULL)
save(MH.alpha.mu.sigma.20.unknown.sn.80.Prawn,file="MH.alpha.mu.sigma.20.unknown.sn.80.Prawn.Rdata")


#=====================================================#
# Results from MH models: Perfect and imperfect Sensitivity
# Targeted Beta and Leakage Parameters
#=====================================================#

load("MH.alpha.mu.sigma.20.Perfect.Prawn.Rdata")            # MH without cut-off k=0
load("MH.alpha.mu.sigma.20.Perfect.Prawn.TBB.0.Rdata")       # MH with cut-off k=0
load("MH.alpha.mu.sigma.20.known.sn.80.Prawn.Rdata")        # MH without cut-off k=0
load("MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0.Rdata") # MH with cut-off k=0 . Used later for Final Tables
load("MH.alpha.mu.sigma.20.unknown.sn.80.Prawn.Rdata")      # MH without cut-off k=0


# Perfect case ----------------------#

summarise_Mean_CI <- function(x) {
  c(mean = mean(x, na.rm = TRUE),
    quantile(x, probs = c(0.025, 0.975), na.rm = TRUE))
}

E.L.P.L.Draws.perfect.case <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                               alpha=MH.alpha.mu.sigma.20.Perfect.Prawn$target.parameters$alpha_sample,
                                               beta=MH.alpha.mu.sigma.20.Perfect.Prawn$target.parameters$beta_sample,
                                               sensitivity=1.0,
                                               cutoff = 0)
summ.param.perfect.case <- summary_mcmc(E.L.P.L.Draws.perfect.case$alpha,
                                        E.L.P.L.Draws.perfect.case$beta,
                                        E.L.P.L.Draws.perfect.case$mu,
                                        E.L.P.L.Draws.perfect.case$E.L,
                                        E.L.P.L.Draws.perfect.case$P.L,
                                        varnames = c("alpha","beta","mu","E.L","P.L"))
round(summ.param.perfect.case$summary_mcmc[,c("Mean","2.5%","97.5%")],5)
# Mean     2.5%    97.5%
# alpha  0.00500  0.00320  0.00733
# beta   6.37092  3.09414 10.52888
# mu     0.00083  0.00051  0.00134
# E.L   22.04209 14.56607 31.33066
# P.L    0.04064  0.02639  0.05865

E.L.P.L.Draws.perfect.case.exact <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                         alpha=unlist((MH.alpha.mu.sigma.20.Perfect.Prawn$target.parameters$alpha_sample)),
                                         beta=unlist((MH.alpha.mu.sigma.20.Perfect.Prawn$target.parameters$beta_sample)),
                                         sensitivity=1.0,
                                         cutoff = 0)
summ.param.perfect.case <- summary_mcmc(E.L.P.L.Draws.perfect.case.exact$alpha,
                                        E.L.P.L.Draws.perfect.case.exact$beta,
                                        E.L.P.L.Draws.perfect.case.exact$mu,
                                        E.L.P.L.Draws.perfect.case.exact$E.L,
                                        E.L.P.L.Draws.perfect.case.exact$P.L,
                                        varnames = c("alpha","beta","mu","E_Li","P_Li"))
round(summ.param.perfect.case$summary_mcmc[,c("Mean","2.5%","97.5%")],5)
# Mean     2.5%    97.5%
# alpha  0.00500  0.00320  0.00733
# beta   6.37092  3.09414 10.52888
# mu     0.00083  0.00051  0.00134
# E_Li  22.04209 14.56607 31.33066
# P_Li   0.04064  0.02639  0.05865


MH.Estimate.Perfect.Case <- summ.param.perfect.case$summary_mcmc[,c("Mean","2.5%","97.5%")]
MH.Estimate.Perfect.Case
#             Mean         2.5%        97.5%
# alpha  0.005002305 3.196561e-03  0.007330628
# beta   6.370922528 3.094140e+00 10.528878986
# mu     0.000828876 5.096304e-04  0.001340304
# E.L   22.042087813 1.456607e+01 31.330661332
# P.L    0.040642436 2.639281e-02  0.058646555

# Using Results with cutoff=0 -------------#

E.L.P.L.Draws.perfect.case.exact <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                                     alpha=unlist((MH.alpha.mu.sigma.20.Perfect.Prawn.TBB.0$target.parameters$alpha_sample)),
                                                     beta=unlist((MH.alpha.mu.sigma.20.Perfect.Prawn.TBB.0$target.parameters$beta_sample)),
                                                     sensitivity=1.0,
                                                     cutoff = 0)
summarise_Mean_CI(E.L.P.L.Draws.perfect.case.exact$E.L)
summarise_Mean_CI(E.L.P.L.Draws.perfect.case.exact$P.L)
summ.param.perfect.case <- summary_mcmc(E.L.P.L.Draws.perfect.case.exact$alpha,
                                        E.L.P.L.Draws.perfect.case.exact$beta,
                                        E.L.P.L.Draws.perfect.case.exact$mu,
                                        E.L.P.L.Draws.perfect.case.exact$E.L,
                                        E.L.P.L.Draws.perfect.case.exact$P.L,
                                        varnames = c("alpha","beta","mu","E.L","P.L"))
round(summ.param.perfect.case$summary_mcmc[,c("Mean","2.5%","97.5%")],5)
# Mean     2.5%    97.5%
#   alpha  0.00506  0.00322  0.00743
# beta   6.65777  3.17491 11.20649
# mu     0.00081  0.00049  0.00132
# E_Li  22.19664 14.64824 31.54576
# P_Li   0.04107  0.02656  0.05932

MH.Estimate.Perfect.Case <- summ.param.perfect.case$summary_mcmc[,c("Mean","2.5%","97.5%")]


# Imperfect case: Known fixed sensitivity : 0.80 ----------------------#

# Using model without cutoff=0 ===========#
E.L.P.L.Draws.imperfect.known.sn <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                             alpha=MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$alpha_sample,
                                             beta=MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$beta_sample,
                                             sensitivity=MH.alpha.mu.sigma.20.known.sn.80.Prawn$target.parameters$sensitivity_sample,
                                             cutoff = 0)
summ.param.imperfect.known.sn <- summary_mcmc(E.L.P.L.Draws.imperfect.known.sn$alpha,
                                            E.L.P.L.Draws.imperfect.known.sn$beta,
                                            E.L.P.L.Draws.imperfect.known.sn$mu,
                                            E.L.P.L.Draws.imperfect.known.sn$E.L,
                                            E.L.P.L.Draws.imperfect.known.sn$P.L,
                                            varnames = c("alpha","beta","mu","E.L","P.L"))
round(summ.param.imperfect.known.sn$summary_mcmc[,c("Mean","2.5%","97.5%")],5)
# Mean     2.5%    97.5%
#   alpha  0.00459  0.00289  0.00681
# beta   4.16470  1.91664  7.33443
# mu     0.00118  0.00070  0.00196
# E.L   25.71597 16.75710 36.94907
# P.L    0.03645  0.02220  0.05436

# Using Model with cutoff=0 ========#

E.L.P.L.Draws.imperfect.known.sn.k.0 <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                               alpha=MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0$target.parameters$alpha_sample,
                                               beta=MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0$target.parameters$beta_sample,
                                               sensitivity=MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0$target.parameters$sensitivity_sample,
                                               cutoff = 0)
summ.param.imperfect.known.sn.k.0 <- summary_mcmc(E.L.P.L.Draws.imperfect.known.sn.k.0$alpha,
                                                  E.L.P.L.Draws.imperfect.known.sn.k.0$beta,
                                                  E.L.P.L.Draws.imperfect.known.sn.k.0$mu,
                                                  E.L.P.L.Draws.imperfect.known.sn.k.0$E.L,
                                                  E.L.P.L.Draws.imperfect.known.sn.k.0$P.L,
                                              varnames = c("alpha","beta","mu","E.L","P.L"))
round(summ.param.imperfect.known.sn.k.0$summary_mcmc[,c("Mean","2.5%","97.5%")],5)
# Mean     2.5%    97.5%
#   alpha  0.00464  0.00291  0.00692
# beta   4.34128  1.96651  7.82616
# mu     0.00115  0.00066  0.00192
# E.L   25.90651 16.78505 37.30170
# P.L    0.03690  0.02242  0.05522

# We use "MH.alpha.mu.sigma.20.Imperfect.80.Prawn.TBB.0" for consistency in Tables



# Imperfect case: unKnown fixed sensitivity : U(0.75, 0.85) ----------------------#
# Here we used without cutoff=0.
# We can update this later if require

E.L.P.L.Draws.imperfect.unknown.sn <- estimate_leakage(B=8000,b=13,M=40,m=5,
                                               alpha=MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$alpha_sample,
                                               beta=MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$beta_sample,
                                               sensitivity=MH.alpha.mu.sigma.20.unknown.sn.80.Prawn$target.parameters$sensitivity_sample,
                                               cutoff=0.0)
MH.Estimate.imperfect.unknown.sn <- summary_mcmc(E.L.P.L.Draws.imperfect.unknown.sn$alpha,
                                              E.L.P.L.Draws.imperfect.unknown.sn$beta,
                                              E.L.P.L.Draws.imperfect.unknown.sn$mu,
                                              E.L.P.L.Draws.imperfect.unknown.sn$E.L,
                                              E.L.P.L.Draws.imperfect.unknown.sn$P.L,
                                              varnames = c("alpha","beta","mu","E.L","P.L"))
MH.Estimate.imperfect.unknown.sn$summary_mcmc[,c("Mean","2.5%","97.5%")]


# Combined Results: MLE, BRMS & MH ---------------

alpha_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"alpha"],pllf.CI.alpha.sn.1$mle.ll,pllf.CI.alpha.sn.1$mle.ul)
alpha_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"alpha"],pllf.CI.alpha.sn.0.80$mle.ll,pllf.CI.alpha.sn.0.80$mle.ul)

beta_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"beta"],pllf.CI.beta.sn.1$mle.ll,pllf.CI.beta.sn.1$mle.ul)
beta_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"beta"],pllf.CI.beta.sn.0.80$mle.ll,pllf.CI.beta.sn.0.80$mle.ul)

mu_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"mu"],pllf.CI.mu.sn.1$mle.ll,pllf.CI.mu.sn.1$mle.ul)
mu_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"mu"],pllf.CI.mu.sn.0.80$mle.ll,pllf.CI.mu.sn.0.80$mle.ul)

E.L_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"E.L"],pllf.CI.E.L.sn.1$mle.ll,pllf.CI.E.L.sn.1$mle.ul)
E.L_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"E.L"],pllf.CI.E.L.sn.0.80$mle.ll,pllf.CI.E.L.sn.0.80$mle.ul)

P.L_estimate_SN.1_MLE <- c(MLE.perfect.case[1,"P.L"],pllf.CI.P.L.sn.1$mle.ll,pllf.CI.P.L.sn.1$mle.ul)
P.L_estimate_SN.0.80_MLE <- c(MLE.imperfect.case.sn.80[1,"P.L"],pllf.CI.P.L.sn.0.80$mle.ll,pllf.CI.P.L.sn.0.80$mle.ul)


alpha_estimate_BRMS <- alpha_estimate_prawn
beta_start_estimate_BRMS <- beta_estimate_prawn   # beta_star=beta/(m*sensitivity)

beta_estimate_SN.1_BRMS <-beta_estimate_prawn_Delta_1.0 # beta=beta_star*m*sensitivity
beta_estimate_SN.0.80_BRMS <-beta_estimate_prawn_Delta_0.80 # beta=beta_star*m*sensitivity

mu_estimate_SN.1_BRMS <-mu_estimate_prawn_Delta_1.0
mu_estimate_SN.0.80_BRMS <-mu_estimate_prawn_Delta_0.80

E.L_estimate_SN.1_BRMS <-E.L_estimate_prawn_Delta_1
E.L_estimate_SN.0.80_BRMS <-E.L_estimate_prawn_Delta_0.80

P.L_estimate_SN.1_BRMS <-P.L_estimate_prawn_Delta_1
P.L_estimate_SN.0.80_BRMS <-P.L_estimate_prawn_Delta_0.80


SN.1.0.MLE <- c(alpha_estimate_SN.1_MLE,beta_estimate_SN.1_MLE,mu_estimate_SN.1_MLE,E.L_estimate_SN.1_MLE,P.L_estimate_SN.1_MLE)
SN.1.0.MH <- c(summ.param.perfect.case$summary_mcmc["alpha",c("Mean","2.5%","97.5%")],
               summ.param.perfect.case$summary_mcmc["beta",c("Mean","2.5%","97.5%")],
               summ.param.perfect.case$summary_mcmc["mu",c("Mean","2.5%","97.5%")],
               summ.param.perfect.case$summary_mcmc["E.L",c("Mean","2.5%","97.5%")],
               summ.param.perfect.case$summary_mcmc["P.L",c("Mean","2.5%","97.5%")])
SN.1.0.BRMS <- c(alpha_estimate_BRMS,beta_estimate_SN.1_BRMS,mu_estimate_SN.1_BRMS,E.L_estimate_SN.1_BRMS,P.L_estimate_SN.1_BRMS)

SN.0.80.MLE <- c(alpha_estimate_SN.0.80_MLE,beta_estimate_SN.0.80_MLE,mu_estimate_SN.0.80_MLE,E.L_estimate_SN.0.80_MLE,P.L_estimate_SN.0.80_MLE)
SN.0.80.MH <- c(MH.Estimate.imperfect.known.sn.80$summary_mcmc["alpha",c("Mean","2.5%","97.5%")],
                MH.Estimate.imperfect.known.sn.80$summary_mcmc["beta",c("Mean","2.5%","97.5%")],
                MH.Estimate.imperfect.known.sn.80$summary_mcmc["mu",c("Mean","2.5%","97.5%")],
                MH.Estimate.imperfect.known.sn.80$summary_mcmc["E.L",c("Mean","2.5%","97.5%")],
                MH.Estimate.imperfect.known.sn.80$summary_mcmc["P.L",c("Mean","2.5%","97.5%")])
SN.0.80.BRMS <- c(alpha_estimate_BRMS,beta_estimate_SN.0.80_BRMS,mu_estimate_SN.0.80_BRMS,E.L_estimate_SN.0.80_BRMS,P.L_estimate_SN.0.80_BRMS)

Unknown.SN.0.80.MH <- c(MH.Estimate.imperfect.unknown.sn$summary_mcmc["alpha",c("Mean","2.5%","97.5%")],
                        MH.Estimate.imperfect.unknown.sn$summary_mcmc["beta",c("Mean","2.5%","97.5%")],
                        MH.Estimate.imperfect.unknown.sn$summary_mcmc["mu",c("Mean","2.5%","97.5%")],
                        MH.Estimate.imperfect.unknown.sn$summary_mcmc["E.L",c("Mean","2.5%","97.5%")],
                        MH.Estimate.imperfect.unknown.sn$summary_mcmc["P.L",c("Mean","2.5%","97.5%")])

Results.Exact.Approx <- rbind(SN.1.0.MLE,
                              SN.1.0.MH,
                              SN.1.0.BRMS,
                              SN.0.80.MLE,
                              SN.0.80.MH,
                              SN.0.80.BRMS,
                              Unknown.SN.0.80.MH)
colnames(Results.Exact.Approx) <- c("alpha","alpha.ll","alpha.ul",
                                    "beta","beta.ll","beta.ul",
                                    "mu","mu.ll","mu.ul",
                                    "E.L","E.L.ll","E.L.ul",
                                    "P.L","P.L.ll","P.L.ul")
rownames(Results.Exact.Approx) <- c(rep("SN 1.0",3),rep("SN 0.80",3),"SN ~ U(0.75,0.85)")

Results.Exact.Approx <- as.data.frame(Results.Exact.Approx)
Results.Exact.Approx$Scenario <- c(rep("SN 1.0",3),rep("SN 0.80",3),"SN ~ U(0.75,0.85)")
Results.Exact.Approx$Method <- c(rep(c("ML","MH","BRMS"),2),rep("MH",1))
Results.Exact.Approx$Delta <- c(rep("1.0",3),rep("0.80",3),"U(0.75,0.85)")
Results.Exact.Approx$Lambda <- c(rep("1.0",3),rep("1.0",3),"1.0")

Results.Exact.Approx <- Results.Exact.Approx[,c("Scenario","Delta","Lambda","Method",
                                                "alpha","alpha.ll","alpha.ul",
                                                "beta","beta.ll","beta.ul",
                                                "mu","mu.ll","mu.ul",
                                                "E.L","E.L.ll","E.L.ul",
                                                "P.L","P.L.ll","P.L.ul")]

Results.Exact.Approx$alpha <- Results.Exact.Approx$alpha*10^3
Results.Exact.Approx$alpha.ll <- Results.Exact.Approx$alpha.ll*10^3
Results.Exact.Approx$alpha.ul <- Results.Exact.Approx$alpha.ul*10^3

Results.Exact.Approx$mu <- Results.Exact.Approx$mu*10^3
Results.Exact.Approx$mu.ll <- Results.Exact.Approx$mu.ll*10^3
Results.Exact.Approx$mu.ul <- Results.Exact.Approx$mu.ul*10^3


Results.Exact.Approx$P.L <- Results.Exact.Approx$P.L*10^2
Results.Exact.Approx$P.L.ll <- Results.Exact.Approx$P.L.ll*10^2
Results.Exact.Approx$P.L.ul <- Results.Exact.Approx$P.L.ul*10^2

rownames(Results.Exact.Approx) <- NULL

write.csv(Results.Exact.Approx,file="ML.MH.BRMS.Results.Exact.Approx.JASA.csv")
save(Results.Exact.Approx,file="ML.MH.BRMS.Results.Exact.Approx.JASA.Rdata")

# Create Table 3 : Main Manuscript ----------------

load("ML.MH.BRMS.Results.Exact.Approx.JASA.Rdata")
library(xtable)
# Table 2
print(xtable((Results.Exact.Approx[,c("Delta","Lambda","Method","alpha","alpha.ll","alpha.ul","beta","beta.ll","beta.ul","mu","mu.ll","mu.ul")]), digits = 2), include.rownames=FALSE)

# Table 3
print(xtable((Results.Exact.Approx[,c("Delta","Lambda","Method","E.L","E.L.ll","E.L.ul","P.L","P.L.ll","P.L.ul")]), digits = 4), include.rownames=FALSE)


# Combine Table 2 and 3

str(Results.Exact.Approx)

library(dplyr)
library(tidyr)

# Reshape and format the desired parameters
param_table <- Results.Exact.Approx %>%
  mutate(
    alpha = sprintf("%.2f (%.2f, %.2f)", alpha, alpha.ll, alpha.ul),
    beta  = sprintf("%.2f (%.2f, %.2f)", beta, beta.ll, beta.ul),
    mu    = sprintf("%.2f (%.2f, %.2f)", mu, mu.ll, mu.ul),
    E.L   = sprintf("%.2f (%.2f, %.2f)", E.L, E.L.ll, E.L.ul),
    P.L   = sprintf("%.2f (%.2f, %.2f)", P.L, P.L.ll, P.L.ul)
  ) %>%
  select(Scenario, Method, alpha, beta, mu, E.L, P.L) %>%
  pivot_longer(cols = alpha:P.L, names_to = "Parameter", values_to = "Estimate") %>%
  pivot_wider(names_from = Method, values_from = Estimate) %>%
  arrange(Scenario, factor(Parameter, levels = c("alpha", "beta", "mu", "E.L", "P.L")))

print(xtable(param_table, digits = 2), include.rownames=FALSE)


