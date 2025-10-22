#-------------------------------------------------------------@
# Compiling Results for Figure 1: JASA PAper
# Using Prawn data 
# Exact BB Model: Frequentist method using MLE and PLLF 
# Using Generailzed Function
# Target Parameters: PLLF of Sensitivity and Specificity
# Relevant Supplementary Table will be from here
#-------------------------------------------------------------@

source("R Functions JASA Paper Updated September 2025.R")

# Data -----

df.prawn <- read.csv("deidentified_frozen_seafood.csv")
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


# PLLF using Final Generalized Function ----------------

# Specificity under imperfect sensitivity

store_specificity_values_prawn_imperfect <- function(specificity) {
  # Run the optimization
  optim_result <- optim(c(0, 0), loglikfn.JABES.sp.sn.TrBB, 
                        ty = summ.prawn.data$ty, 
                        freq = summ.prawn.data$freq, 
                        b = 13, 
                        m = 5,
                        M = 40,
                        theta = Inf, 
                        R = 1e4, 
                        # R = 3e4,
                        sensitivity = 0.8, 
                        specificity = specificity, 
                        cutoff = 0,
                        deviance = FALSE, 
                        control = list(reltol = 1e-12, fnscale = -1), 
                        hessian = FALSE)
  
  # Store the parameter estimates
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$cutoff <- 0
  optim_result$mu.TrBB <- optim_result$alpha / (optim_result$alpha + optim_result$beta)
  #  optim_result$mu.overall <- (1 - optim_result$pi) * optim_result$mu.BB
  optim_result$sensitivity <- 0.8
  optim_result$specificity <- specificity
  
  # Return the relevant values as a data frame
  return(data.frame(cutoff = optim_result$cutoff, 
                    alpha = optim_result$alpha, 
                    beta = optim_result$beta, 
                    mu.BB = optim_result$mu.TrBB, 
                    sensitivity = optim_result$sensitivity, 
                    specificity = optim_result$specificity, 
                    #  mu.overall = optim_result$mu.overall, 
                    value = optim_result$value,
                    convergence=optim_result$convergence
  ))
}

specificity <- seq(0.9989,1,by=0.000001)

# Use lapply to apply the function to each pi value and store results as a list of dataframes
PLLF_Specificity_imperfect_Sensitivity_Prawn <- lapply(specificity, store_specificity_values_prawn_imperfect)

# Combine the list of dataframes into a single dataframe
PLLF_Specificity_imperfect_Sensitivity_Prawn <- do.call(rbind, PLLF_Specificity_imperfect_Sensitivity_Prawn)

save(PLLF_Specificity_imperfect_Sensitivity_Prawn,file="PLLF_Specificity_imperfect_Sensitivity_Prawn_Re4.Rdata")


# Sensitivity under perfect Specificity  

store_sensitivity_values_prawn_perfect <- function(sensitivity) {
  # Run the optimization
  optim_result <- optim(c(0, 0), loglikfn.JABES.sp.sn.TrBB, 
                        ty = summ.prawn.data$ty, 
                        freq = summ.prawn.data$freq, 
                        b = 13, 
                        m = 5,
                        M = 40,
                        theta = Inf, 
                        R = 1e4, 
                        # R = 3e4,
                        sensitivity = sensitivity, 
                        specificity = 1.0, 
                        cutoff = 0,
                        deviance = FALSE, 
                        control = list(reltol = 1e-12, fnscale = -1), 
                        hessian = FALSE)
  
  # Store the parameter estimates
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$cutoff <- 0
  optim_result$mu.TrBB <- optim_result$alpha / (optim_result$alpha + optim_result$beta)
  #  optim_result$mu.overall <- (1 - optim_result$pi) * optim_result$mu.BB
  optim_result$sensitivity <- sensitivity
  optim_result$specificity <- 1.0
  
  # Return the relevant values as a data frame
  return(data.frame(cutoff = optim_result$cutoff, 
                    alpha = optim_result$alpha, 
                    beta = optim_result$beta, 
                    mu.BB = optim_result$mu.TrBB, 
                    sensitivity = optim_result$sensitivity, 
                    specificity = optim_result$specificity, 
                    #  mu.overall = optim_result$mu.overall, 
                    value = optim_result$value,
                    convergence=optim_result$convergence
  ))
}


sensitivity  <- seq(0.45,1,by=0.0001)

# Use lapply to apply the function to each pi value and store results as a list of dataframes
PLLF_Sensitivity_perfect_Specificity_Prawn <- lapply(sensitivity, store_sensitivity_values_prawn_perfect)

# Combine the list of dataframes into a single dataframe
PLLF_Sensitivity_perfect_Specificity_Prawn <- do.call(rbind, PLLF_Sensitivity_perfect_Specificity_Prawn)

save(PLLF_Sensitivity_perfect_Specificity_Prawn,file="PLLF_Sensitivity_perfect_Specificity_Prawn_Re4.Rdata")



# Profile likelihood CI: specificity under perfect sensitivity --------------
# Specificity vary given Sensitivty=1

specificity <- seq(0.99,1,by=0.00001)
sensitivity  <- 1
df.pllf.sp <- expand.grid(sensitivity,specificity)
dim(df.pllf.sp)

colnames(df.pllf.sp) <- c("sensitivity","specificity")
df.pllf.sp$logL <- rep(NA,dim(df.pllf.sp)[1])
df.pllf.sp$alpha <- rep(NA,dim(df.pllf.sp)[1])
df.pllf.sp$beta <- rep(NA,dim(df.pllf.sp)[1])
df.pllf.sp$mu <- rep(NA,dim(df.pllf.sp)[1])
df.pllf.sp$converge <- rep(NA,dim(df.pllf.sp)[1])
df.pllf.sp$P.L <- rep(NA,dim(df.pllf.sp)[1])
df.pllf.sp$E.L <- rep(NA,dim(df.pllf.sp)[1])



for (i in 1:dim(df.pllf.sp)[1]){
  
  skip_to_next <- FALSE
  
  tryCatch(
    
    {o.theta.inf.prawn.sp <- optim(c(0,0),loglikfn.JABES.sp.sn,ty=c(summ.prawn.data$ty),freq=c(summ.prawn.data$freq),b=13,m=5,
                                   theta=Inf,R=1e4,sensitivity=df.pllf.sp$sensitivity[i],specificity=df.pllf.sp$specificity[i],
                                   deviance=FALSE,control=list(reltol=1e-12,fnscale=-1),hessian=TRUE)
    
    temp <- loglikfn.JABES.sp.sn(o.theta.inf.prawn.sp$par,ty=summ.prawn.data$ty,freq=summ.prawn.data$freq,b=13,B=8000,
                                 theta=Inf,m=5,M=40,leakage=TRUE,R=1e4,sensitivity=df.pllf.sp$sensitivity[i],specificity=df.pllf.sp$specificity[i])
    o.theta.inf.prawn.sp$P.L <- temp["P.L"]
    o.theta.inf.prawn.sp$E.L <- temp["E.L"]
    
    
    df.pllf.sp$logL[i] <- o.theta.inf.prawn.sp$value
    df.pllf.sp$mu[i] <- exp(o.theta.inf.prawn.sp$par[1])/sum(exp(o.theta.inf.prawn.sp$par))
    df.pllf.sp$alpha[i] <- exp(o.theta.inf.prawn.sp$par[1])
    df.pllf.sp$beta[i] <- exp(o.theta.inf.prawn.sp$par[2])
    df.pllf.sp$converge[i] <- o.theta.inf.prawn.sp$convergence
    df.pllf.sp$P.L[i] <- o.theta.inf.prawn.sp$P.L
    df.pllf.sp$E.L[i] <- o.theta.inf.prawn.sp$E.L
    },
    error = function(e) { skip_to_next <<- TRUE})
  
  if(skip_to_next) { next }
  
  
}


df.pllf.sp.sn.1.prawn <- df.pllf.sp

save(df.pllf.sp.sn.1.prawn,file="df.pllf.sp.sn.1.prawn.rdata")





# Figure 1: Manuscript -------------------------------


load("df.pllf.sp.sn.1.prawn.rdata") # Using old function
load("df.pllf.sn.sp.1.prawn.rdata") # Using old function

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



# Figure 1 ------

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


# Figure 1 : Re-scaled ------

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



# Supplementary Material Analysis : Run PLLF --------------


# PLLF for Specificity: Without restriction on sensitivity -------------------#
# Using Previous function: loglikfn.JABES.sp.sn

o.theta.inf.prawn.sp.sn <- optim(c(0,0,0),loglikfn.JABES.sp.sn,ty=summ.prawn.data$ty,freq=summ.prawn.data$freq,b=13,
                                 theta=Inf,m=5,M=40,specificity=1,deviance=FALSE,control=list(factr=1e-12,fnscale=-1,trace=T),
                                 R=1e4,hessian=FALSE,method = "L-BFGS-B",lower = c(log(0.0001),log(0.0001),log(0.50)),upper = c(Inf,Inf,0))

store_specificity_values_prawn_imperfect <- function(specificity) {
  # Run the optimization
  optim_result <- optim(c(0, 0, 0), loglikfn.JABES.sp.sn, 
                        ty = summ.prawn.data$ty, 
                        freq = summ.prawn.data$freq, 
                        b = 13, 
                        m = 5,
                        M = 40,
                        theta = Inf, 
                        R = 1e4, 
                        # R = 3e4,
                        # sensitivity = 0.8, 
                        specificity = specificity, 
                        deviance = FALSE,
                        hessian = FALSE,
                        control=list(factr=1e-12,fnscale=-1,trace=T),
                        method = "L-BFGS-B",
                        lower = c(log(0.0001),log(0.0001),log(0.10)),
                        upper = c(Inf,Inf,0)
                        )
  
  # Store the parameter estimates
  optim_result$alpha <- exp(optim_result$par[1])
  optim_result$beta <- exp(optim_result$par[2])
  optim_result$sensitivity <- exp(optim_result$par[3])
  
  optim_result$mu <- optim_result$alpha / (optim_result$alpha + optim_result$beta)
  #  optim_result$mu.overall <- (1 - optim_result$pi) * optim_result$mu.BB
  optim_result$specificity <- specificity
  
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

specificity <- seq(0.9989,1,by=0.0000005)
# specificity <- seq(0.9989,1,by=0.0001)

# Use lapply to apply the function to each pi value and store results as a list of dataframes

PLLF_Specificity_Free_Parameter_Prawn <- lapply(specificity, function(spcs) {
  tryCatch({
    store_specificity_values_prawn_imperfect(spcs)
  }, error = function(e) {
    message("Error for specificity = ", spcs, ": ", e$message)
    return(NA)  # or NULL, depending on your downstream needs
  })
})




# Combine the list of dataframes into a single dataframe
PLLF_Specificity_Free_Parameter_Prawn <- do.call(rbind, PLLF_Specificity_Free_Parameter_Prawn)

# save(PLLF_Specificity_Free_Parameter_Prawn,file="PLLF_Specificity_Free_Parameter_Prawn_Re4.Rdata")




# PLLF for Sensitivity: Without restriction on Specificity -------------------#
# Using Previous function: loglikfn.JABES.sp.sn

o.theta.inf.prawn.sp.sn <- optim(c(0,0,0),loglikfn.JABES.sp.sn,ty=summ.prawn.data$ty,freq=summ.prawn.data$freq,b=13,
                                 theta=Inf,m=5,M=40,sensitivity=1,deviance=FALSE,control=list(factr=1e-12,fnscale=-1,trace=T),
                                 R=1e4,hessian=FALSE,method = "L-BFGS-B",lower = c(log(0.0001),log(0.0001),log(0.9999)),upper = c(Inf,Inf,0))

store_sensitivity_values_prawn_imperfect <- function(sensitivity) {
  # Run the optimization
  optim_result <- optim(c(0, 0, 0), loglikfn.JABES.sp.sn, 
                        ty = summ.prawn.data$ty, 
                        freq = summ.prawn.data$freq, 
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
  #  optim_result$mu.overall <- (1 - optim_result$pi) * optim_result$mu.BB
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
# sensitivity  <- seq(0.45,1,by=0.001)

# Use lapply to apply the function to each pi value and store results as a list of dataframes

# PLLF_Sensitivity_Free_Parameter_Prawn <- lapply(sensitivity, store_sensitivity_values_prawn_imperfect)

PLLF_Sensitivity_Free_Parameter_Prawn <- lapply(sensitivity, function(sens) {
  tryCatch({
    store_sensitivity_values_prawn_imperfect(sens)
  }, error = function(e) {
    message("Error for sensitivity = ", sens, ": ", e$message)
    return(NA)  # or NULL, depending on your downstream needs
  })
})



# Combine the list of dataframes into a single dataframe
PLLF_Sensitivity_Free_Parameter_Prawn <- do.call(rbind, PLLF_Sensitivity_Free_Parameter_Prawn)

save(PLLF_Sensitivity_Free_Parameter_Prawn,file="PLLF_Sensitivity_Free_Parameter_Prawn_Re4_small_scale.Rdata")

PLLF_df <- as.data.frame(t(PLLF_Sensitivity_Free_Parameter_Prawn))
colnames(PLLF_df) <- rownames(PLLF_Sensitivity_Free_Parameter_Prawn)

# save(PLLF_Sensitivity_Free_Parameter_Prawn,file="PLLF_Sensitivity_Free_Parameter_Prawn_Re4.Rdata")


# Create Plots -------------------------------------#

load("PLLF_Specificity_Free_Parameter_Prawn_Re4.Rdata")

PLLF_Specificity_Free_Parameter_Prawn <- PLLF_Specificity_Free_Parameter_Prawn[!is.na(PLLF_Specificity_Free_Parameter_Prawn$convergence) & PLLF_Specificity_Free_Parameter_Prawn$convergence==0,]
PLLF_Specificity_Free_Parameter_Prawn$logL <- PLLF_Specificity_Free_Parameter_Prawn$value

pllf.CI.estimate.specificity.95 <- pllf.CI.estimate(PLLF_Specificity_Free_Parameter_Prawn$value,parameter = PLLF_Specificity_Free_Parameter_Prawn$specificity,level = 0.95)
pllf.CI.estimate.specificity.90 <- pllf.CI.estimate(PLLF_Specificity_Free_Parameter_Prawn$value,parameter = PLLF_Specificity_Free_Parameter_Prawn$specificity,level = 0.90)

plot(PLLF_Specificity_Free_Parameter_Prawn$specificity,PLLF_Specificity_Free_Parameter_Prawn$value,typ="l")
plot(PLLF_Specificity_Free_Parameter_Prawn$sensitivity,PLLF_Specificity_Free_Parameter_Prawn$value,typ="l")





# load("PLLF_Sensitivity_Free_Parameter_Prawn_Re4_small_scale.Rdata")
load("PLLF_Sensitivity_Free_Parameter_Prawn_Re4.Rdata")

PLLF_Sensitivity_Free_Parameter_Prawn <- PLLF_Sensitivity_Free_Parameter_Prawn[!is.na(PLLF_Sensitivity_Free_Parameter_Prawn$convergence) & PLLF_Sensitivity_Free_Parameter_Prawn$convergence==0,]

o.theta.inf.prawn.sn.sp.1 <- optim(c(0,0),loglikfn.JABES.sp.sn,ty=summ.prawn.data$ty,freq=summ.prawn.data$freq,b=13,
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



# Plot: On original scale --------#

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


# Plot by Rescaling: LogL - LogL at MLE (in Supplementary File) ------------------# 




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


# Examine the Sensitivity under imperfect specificity: Not in Supplementary File --------------

load("PLLF_Sensitivity_Free_Parameter_Prawn_Re4.Rdata")

PLLF_Sensitivity_Free_Parameter_Prawn <- PLLF_Sensitivity_Free_Parameter_Prawn[!is.na(PLLF_Sensitivity_Free_Parameter_Prawn$convergence) & PLLF_Sensitivity_Free_Parameter_Prawn$convergence==0,]
PLLF_Sensitivity_Free_Parameter_Prawn$logL <- PLLF_Sensitivity_Free_Parameter_Prawn$value

write.csv(PLLF_Sensitivity_Free_Parameter_Prawn,file="PLLF_Sensitivity_Free_Parameter_Prawn.csv")

PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc <- PLLF_Sensitivity_Free_Parameter_Prawn[PLLF_Sensitivity_Free_Parameter_Prawn$specificity==1,]
PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc <- PLLF_Sensitivity_Free_Parameter_Prawn[PLLF_Sensitivity_Free_Parameter_Prawn$specificity<1,]

pllf.CI.estimate.sensitivity.95 <- pllf.CI.estimate(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$value,parameter = PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$sensitivity,level = 0.95)
pllf.CI.estimate.sensitivity.90 <- pllf.CI.estimate(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$value,parameter = PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$sensitivity,level = 0.90)

png(file="pllf CI sensitivity free parameter.png", width = 12,height = 4,units = "in",res = 300)

par(mfrow=c(1,3), mgp=c(3, 0.5, 0), oma=c(0,0,0,0), mar=c(5, 4, 4, 2))  # Adjusted margins


range_ll <- range(PLLF_Sensitivity_Free_Parameter_Prawn$logL)

vals.2 <- seq(round(range_ll[1],1), round(range_ll[2],2), length.out = 15)
formatted_labels_2 <- formatC(vals.2, format = "f", digits = 2)
vals.1 <- round(seq(min(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity),max(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity), length.out = 15),2)
formatted_labels_1 <- formatC(vals.1, format = "f", digits = 2)

plot(PLLF_Sensitivity_Free_Parameter_Prawn$sensitivity,PLLF_Sensitivity_Free_Parameter_Prawn$logL,typ="p",xlab="Sensitivity",ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1,
     main="logL against sensitivity") # 
axis(side = 2, at = vals.2, labels = formatted_labels_2 ,las = 1,cex.lab=0.6)
axis(side = 1, at = vals.1, labels = formatted_labels_1 , las = 2,cex.lab=0.6)

range_ll <- range(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$logL)

vals.2 <- seq(round(range_ll[1],1), round(range_ll[2],2), length.out = 15)
formatted_labels_2 <- formatC(vals.2, format = "f", digits = 2)
vals.1 <- round(seq(min(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$sensitivity),max(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$sensitivity), length.out = 15),2)
formatted_labels_1 <- formatC(vals.1, format = "f", digits = 2)

plot(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$sensitivity,PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$logL,typ="l",xlab="Sensitivity",ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1,
     main="logL: specificity=1.0") # 
axis(side = 2, at = vals.2, labels = formatted_labels_2 ,las = 1,cex.lab=0.6)
axis(side = 1, at = vals.1, labels = formatted_labels_1 , las = 2,cex.lab=0.6)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$logL,na.rm=T),col=2,lty=3)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$logL,na.rm=T)-3.84/2,col=4)
abline(h=max(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$logL,na.rm=T)-2.71/2,col=5)
abline(v=pllf.CI.estimate.sensitivity.95$mle,col=2,lty=1)
abline(v=pllf.CI.estimate.sensitivity.95$mle.ll,col=4,lty=3)
abline(v=pllf.CI.estimate.sensitivity.95$mle.ul,col=4,lty=3)
abline(v=pllf.CI.estimate.sensitivity.90$mle.ll,col=5,lty=6)
abline(v=pllf.CI.estimate.sensitivity.90$mle.ul,col=5,lty=6)
# legend("bottomright",legend=c(paste("MLE=",mle.sensitivity),paste("LL=",ll.sensitivity),paste("UL=",ul.sensitivity)),bty = "n")


# Identify the x-value corresponding to the maximum y-value
max_x <- pllf.CI.estimate.sensitivity.95$mle
max_y <- max(PLLF_Sensitivity_Free_Parameter_Prawn_perf_spc$logL,na.rm=T)

# Add a vertical line at the max x-value
#abline(v=max_x, col="red", lty=2)

# Display the x-value on the plot
text(x=max_x-0.01, y=(max_y-1.5), labels=round(max_x, 5), pos=3, col=2, cex=0.7, srt = 90)
text(x=max_x+0.25, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.95$mle.ul, 5), pos=3, col=6, cex=0.5, srt = 90)
text(x=max_x-0.27, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.95$mle.ll, 5), pos=3, col=4, cex=0.5, srt = 90)
text(x=max_x-0.235, y=(max_y-1.5), labels=round(pllf.CI.estimate.sensitivity.90$mle.ll, 5), pos=3, col=5, cex=0.5, srt = 90)

legend("bottomright",legend=c("0.95","0.90"),col=c(4,5),lty=c(3,6),bty="n",lwd=c(2,2))

range_ll <- range(PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$logL)

vals.2 <- seq(round(range_ll[1],1), round(range_ll[2],2), length.out = 15)
formatted_labels_2 <- formatC(vals.2, format = "f", digits = 2)
vals.1 <- round(seq(min(PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$sensitivity),max(PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$sensitivity), length.out = 15),5)
formatted_labels_1 <- formatC(vals.1, format = "f", digits = 2)

sensitivity_sorted <- PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$sensitivity[order(PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$specificity)]
specificity_sorted <- PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$specificity[order(PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$specificity)]
logL_sorted <- PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$logL[order(PLLF_Sensitivity_Free_Parameter_Prawn_imperf_spc$specificity)]
plot(sensitivity_sorted,
     logL_sorted,
     typ="p",xlab="Sensitivity",ylab="logL",yaxt="n",xaxt="n",cex.lab=1.2,cex.axis=1,
     main="logL: specificity < 1.0") # 
axis(side = 2, at = vals.2, labels = formatted_labels_2 ,las = 1,cex.lab=0.6)
axis(side = 1, at = vals.1, labels = formatted_labels_1 , las = 2,cex.lab=0.6)

# Add specificity as text labels
text(sensitivity_sorted, logL_sorted, labels = round(specificity_sorted, 6),
     pos = 3, cex = 0.5, col = "blue")  # pos=3 = above the point

dev.off()
