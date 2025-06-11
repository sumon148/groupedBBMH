# Data ------------------------
# Data for test and paper
load("Prawn Data/df.prawn.Rdata")
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

b <- 13
Nbar <- 5
ty <- prawn.data$ty
D <- length(ty)
ty2 <- x
wt <- freq
d <- length(ty2)

# Implementation of R Functions ------------

# Test: rtrbeta --------------#
set.seed(123)
pi.beta <- rbeta(10, shape1 = 2, shape2 = 5)
set.seed(123)
pi.trbeta <-rtrbeta(10, shape1 = 2, shape2 = 5, cutoff = 0)
plot(pi.beta,pi.trbeta)
abline(0,1)

# Test: logit and inverse_logit --------------#
inverse_logit(0.8472979)
logit(0.7)

# Test: MH.Sampler.Chain.TBB  --------------#

G=500;R=1e3;
par <- c(log(0.005), logit(7e-04)) ; initial <- par; k <- 0.0000
MH.alpha.mu.sigma.20.Perfect.Prawn.TBB.0 <- MH.Sampler.Chain.TBB(par=par,cutoff=k,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                                 sigma=0.2,num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=FALSE,specificity=FALSE,
                                                                 sensitivity.range=NULL,specificity.range=NULL)
G=500;R=1e3;
par <- c(log(0.005), logit(7e-04)) ; initial <- par; k <- 0.0025
MH.alpha.mu.sigma.20.Perfect.Prawn.TBB.0025 <- MH.Sampler.Chain.TBB(par=par,cutoff=k,initial=par,alpha=TRUE,beta=FALSE,mu=TRUE,rho=FALSE,G=G,R=R,
                                                                 sigma=0.2,num_chains=4,burnin=R*0.5,thin=5,seed=c(123,456,789,135),trail=b,ty2=ty2,wt=wt,group.size=Nbar,sensitivity=FALSE,specificity=FALSE,
                                                                 sensitivity.range=NULL,specificity.range=NULL)

# Test: E.P.Leakage.trbb  --------------#

alpha_tbb_00 <- unlist(MH.alpha.mu.sigma.20.Perfect.Prawn.TBB.0$target.parameters$alpha_sample)
beta_tbb_00 <- unlist(MH.alpha.mu.sigma.20.Perfect.Prawn.TBB.0$target.parameters$beta_sample)

E.P.trBB.cutoff.00.approx <- E.P.Leakage.trbb(alpha = alpha_tbb_00,beta = beta_tbb_00,cutoff = 0,B = 8000,M = 40,b = 13,m = 5,Delta = 0.70)
hist(E.P.trBB.cutoff.00.approx$alpha)
hist(E.P.trBB.cutoff.00.approx$beta)
hist(E.P.trBB.cutoff.00.approx$E_Li)
hist(E.P.trBB.cutoff.00.approx$P_Li)
