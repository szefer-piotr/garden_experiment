# Random model

#### NOTES -----------
# Simulate random effect model.
# 
# m <- 5       # number of schools
# mu <- 12.56  # global score average
# sdU <- 5.67  # standard deviation for the random effect
# sdW <- 3.25  # standard deviation for the noise
# sU <- sdU^2  # variance for the random effect
# sW <- sdW^2
# # values for the random effect
# U <- rnorm(m, mean=0, sd = sdU)
# # Number of observations in each school
# n <- 100
# N <- m*n
# data <- data.frame(school = rep(1:m, each=n))
# data$U <- U[data$school] # random factor value for each school
# data$val <- rnorm(m*n, mu + data$U, sdW)
# library(lme4)
# library(lmerTest)
# 
# # REML estimator for sigma^2
# MSE <- 1/(m*(n-1)) * sum((data$val - mean(data$val))^2) 
# # REML estimator for sigma_a^2
# MSA <- (1/(m-1))*sum(tapply(data$val, data$school, function(x){sum((x - mean(x))^2)}))
# (1/n)*(MSA - MSE)
# # Estimate parameters set earlier
# m1 <- lmer(val ~ 1 + (1|school), data = data)
# 
# 
# summary(m1)
# MSA
# MSE
# 
# # Predicting random effect
# # Random effect
# estre <- 9.545 
# estvr <- 8.188
# 
# summary(m1)
# 
# i = 5
# a1 = estre/(estre - estvr/n)*(mean(data[data$school==i,"val"]) - mean(data$val))
# data[data$school == i, "U"][1]
# a1


# Mixed model

# ngard <- 6    # number of blocks/gardens
# nplot <- 6
# 
# mu <- 3.75    # global average log[bio]
# fung <- 0.75  # effect size for fungicide
# ins <- -1.5   # effect size for insecticide
# h1 <- 0.5     # effect size herbivory low
# h2 <- 1.0     # effect size herbivory high
# pred <- 0.55  # effect size predator
# 
# sdg <- 0.2791    # standrad deviation for garden (randon effect)
# sd <- 0.4606       # standard deviation for plot
# 
# # Labels for gardens
# # set.seed(16)
# gard = rep(LETTERS[1:ngard], each = nplot)         # garden label
# gardeff = rep( rnorm(ngard, 0, sdg), each = nplot) # garden effect
# ploteff = rnorm(ngard*nplot, 0, sd)                # noise for plots
# 
# # Simulate the explanatory variables
# treat <- c()
# for (g in 1:ngard){treat <- c(treat,
#                               sample(c("c","f","i","h1","h2","p"), nplot))}
# 
# insecticide <- as.numeric(treat == "i")
# fungicide   <- as.numeric(treat == "f")
# herbivory1  <- as.numeric(treat == "h1")
# herbivory2  <- as.numeric(treat == "h2")
# predator    <- as.numeric(treat == "p")
# 
# bio <- mu + fung*fungicide + ins*insecticide + h1*herbivory1 + h2*herbivory2 + pred*predator + gardeff + ploteff
# 
# ds <- as.data.frame(cbind(bio,
#                           gard,
#                           treat))
# 
# ds$bio <- as.numeric(as.character(ds$bio))

# SIMULATIONS -----

library(lme4)
library(lmerTest)
library(ggplot2)
library(lattice)
library(latticeExtra)

generateData <- function(ngard = 6,    # number of blocks/gardens
                         nplot = 6,    
                         mu = 4.40862,     # global average log[bio]
                         fung = -0.19292,  # effect size for fungicide
                         ins = 0.43223,   # effect size for insecticide
                         h1 = -0.07064,     # effect size herbivory low
                         h2 = -0.95572,     # effect size herbivory high
                         pred = -0.47436,  # effect size predator
                         sdg = 0.2791, # standrad dev for the random effect
                         sd = 0.4606   # standard dev for the residuals 
                         ){
  # Labels for gardens
  # set.seed(16)
  
  # Control
  # print(c(paste("randomvar=",sdg^2), 
  #         paste("residualvar=",sd^2),
  #         paste("effsize=",ins)))
  
  gard = rep(LETTERS[1:ngard], each = nplot)         # garden label
  gardeff = rep(rnorm(ngard, 0, sdg^2), each = nplot) # garden effect
  ploteff = rnorm(ngard*nplot, 0, sd^2)                # noise for plots
  
  # Simulate the explanatory variables
  treat <- c()
  for (g in 1:ngard){treat <- c(treat,
                                sample(c("c","f","i","h1","h2","p"), nplot))}
  
  insecticide <- as.numeric(treat == "i")
  fungicide   <- as.numeric(treat == "f")
  herbivory1  <- as.numeric(treat == "h1")
  herbivory2  <- as.numeric(treat == "h2")
  predator    <- as.numeric(treat == "p")
  
  bio <- mu + fung*fungicide + ins*insecticide + h1*herbivory1 + h2*herbivory2 + pred*predator + gardeff + ploteff
  
  ds <- as.data.frame(cbind(bio,
                            gard,
                            treat))
  
  ds$bio <- as.numeric(as.character(ds$bio))
  
  return(ds)
  
}


ds <- generateData(sdg=1, sd=0.5)

plt <- ggplot(ds, aes(x = treat, y =bio))
plt + geom_boxplot() + 
  geom_jitter(aes(color=gard), width = 0.1, size=3)

fitted_model <- lmer(bio~treat+(1|gard), data = ds)
summary(fitted_model)
sqrt(unlist(VarCorr(fitted_model)))

# Decrease the variability between plots and see how the power changes.
# Effect sizes for each treatment... or can this be only for one example treatment, showing how big of a should be to be detected.

pA <- function(rands = 5,...){
  pvals <- c()#store all the pvalues for an effect
  for(rand in 1:rands){
    ds <- generateData(...)
    slds <- summary(lmer(bio~treat+(1|gard), data = ds))
    pvals <- c(pvals,slds$coefficients["treati", 5])
  }
  return(sum(pvals<0.05)/rands)
}

# Create a matrix comparing two different values... for example
# sdgs <-seq(0.5, 2.5, by=0.1)
# pvals <- c()
# for(sdg in sdgs){
#   sim <- pA(rands = 9, sdg = sdg)
#   pvals <- c(pvals, sim)
# }
# 
# plot(sdgs, pvals)
# 
# effs <- seq(0.1, -1.5, by=-0.1)
# pvals <- c()
# for(eff in effs){
#   sim <- pA(rands = 9, ins = eff)
#   pvals <- c(pvals, sim)
# }


resmatgard <- as.matrix(read.table("resmatgard.txt"))

levelplot(resmatgard, xlab = "Effect size % comparing with the mean",
          ylab = "Standard deviation of the random effect") +
  as.layer(contourplot(resmatgard, col = "red"))


# write.table(resmatgard, "resmatgard.txt")


mu = 4.40862
percs <- seq(0.1, 1, by=0.05)
sdvs <- seq(0.1, 1, by = 0.05)

# Random effect standard deviation
resmatgard1 <- matrix(0, nrow = length(percs),
                     ncol = length(sdvs))
colnames(resmatgard1) <- sdvs
rownames(resmatgard1) <- percs
for(perc in percs){
  for(sdv in sdvs){
    sim <- pA(rands = 3, ins = mu*perc, sdg=sdv, sd=0.7, ngard=6)
    resmatgard1[as.character(perc),as.character(sdv)] <- sim
  }
}
levelplot(resmatgard1, xlab = "Effect size % comparing with the mean",
          ylab = "Standard deviation of the random effect") +
  as.layer(contourplot(resmatgard1, col = "red"))
# Statistical power, why is the random effect not important?


# Residual standard deviation
resmatgard2 <- matrix(0, nrow = length(percs),
                     ncol = length(sdvs))
colnames(resmatgard2) <- sdvs
rownames(resmatgard2) <- percs
for(perc in percs){
  for(sdv in sdvs){
    sim <- pA(rands = 3, ins = mu*perc, sdg=0.25, sd=sdv, ngard=6)
    resmatgard2[as.character(perc),as.character(sdv)] <- sim
  }
  
}
levelplot(resmatgard2, xlab = "Effect size % comparing with the mean",
          ylab = "Standard deviation of the random effect") +
  as.layer(contourplot(resmatgard2, col = "red"))


# Residual and random effect
resmatgard3 <- matrix(0, nrow = length(sdvs),
                     ncol = length(sdvs))
colnames(resmatgard3) <- sdvs
rownames(resmatgard3) <- sdvs
for(sdvr in sdvs){
  for(sdvc in sdvs){
    sim <- pA(rands = 3, ins = 0.35*mu, sdg=sdvr, sd=sdvc, ngard=6)
    resmatgard3[as.character(sdvr),as.character(sdvc)] <- sim
  }
}
levelplot(resmatgard3, xlab = "Effect size % comparing with the mean",
          ylab = "Standard deviation of the random effect") +
  as.layer(contourplot(resmatgard3, col = "red"))

## Using simr package ----

# run analyses_wanang_v2.R first
install.packages("simr")
library(simr)
x <- c("TREATWEEVIL125")

fixef(bio_rbl)[x] <- -0.05
powerSim(bio_rbl)

model2 <- extend(bio_rbl, along = "GARDEN", n=20)
powerSim(model2)
pc <- powerCurve(model2)
plot(pc)
