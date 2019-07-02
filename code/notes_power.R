ngard = 9    # number of blocks/gardens
nplot = 6    
mu = 95    # global average log[bio]
fung = 0  # effect size for fungicide
ins = 0   # effect size for insecticide
h1 = 0     # effect size herbivory low
h2 = -5     # effect size herbivory high
pred = 0  # effect size predator
sdg = 7 # standrad dev for the random effect
sd = 7   # standard dev for the residuals 

gard = rep(LETTERS[1:ngard], each = nplot)         # garden label
gardeff = rep(rlnorm(ngard, 0, log(sdg)), each = nplot) # garden effect
ploteff = rlnorm(ngard*nplot, 0, log(sd))                # noise for plots

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

library(ggplot2)
library(lmerTest)
library(simr)

genDat <- function(ngard = 6,    # number of blocks/gardens
                    nplot = 6,
                    mu = 95,   # global average log[bio]
                    fung = 0,  # effect size for fungicide
                    ins = 0,  # effect size for insecticide
                    h1 = 0,     # effect size herbivory low
                    h2 = -5,     # effect size herbivory high
                    pred = 0,  # effect size predator
                    sdg = 0.07332, # standrad dev for the random effect
                    sd = 1){    # standard dev for the residuals
  gard = rep(LETTERS[1:ngard], each = nplot)         # garden label
  gardeff = rep(rlnorm(ngard, 0, sdg), each = nplot) # garden effect
  ploteff = rlnorm(ngard*nplot, 0, sd)                # noise for plots
  # print(paste("Generating data with sd garden:",sdg," and error var", sd))
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

ds <- genDat(ngard=9,ins=-5)

plt <- ggplot(ds, aes(x = treat, y = bio))
plt + geom_boxplot() + 
  geom_jitter(aes(color=gard), width = 0.1, size=3)

# Randomization

randDat <- function(rands=9, pattern=c(1, 0, 0, -15,0,0), sdg=2, sd=2, ngard=6){
  # pattern <- c(int,fung,h1,h2,ins,pred)
  fits <- c()
  for (rand in 1:rands){
    ds <- genDat(ngard=ngard,
                 fung = pattern[2],
                 h1 = pattern[3],
                 h2 = pattern[4],
                 ins=pattern[5],
                 pred = pattern[6],
                 sdg = sdg,
                 sd = sd)
    print(ds)
    fitted_model <- lmer(log(bio)~treat+(1|gard), data = ds)
    sfm <- summary(fitted_model)
    #print(sfm)
    # print(pattern)
    res <- sum((sfm$coefficients[, "Pr(>|t|)"] < 0.05) == (pattern != 0)) == 6
    fits <- c(fits, res)
  }
  return(fits)
}

# pattern <- c(int,fung,h1,h2,ins,pred)
rd <- randDat(rands=1, pattern=c(1, 0, 0, -10,0,0), ngard=9, sdg = 5, sd =5)
sum(rd/999)

p_hat <- sum(rd/mcrands)
# p_hat + c(-1,1)*1.96*sqrt(p_hat*(1-p_hat)/mcrands)

# Bootstraped confidence intervals
bootSE <- sqrt((sum((rd - mean(rd))^2)/(mcrands -1)))
bootVar <- bootSE^2
# p_hat + c(-1,1)*1.96*sqrt(bootVar/mcrands)


fitted_model <- lmer(log(bio)~treat+(1|gard), data = ds)
sfm <- summary(fitted_model)
pattern <- c(1, fung, h1, h2,ins, pred)

# Test for the same pattern
sum(((sfm$coefficients[, "Pr(>|t|)"] < 0.05) == (pattern>0))) == 6



powerSim(fitted_model, VarCorr =)

sigma2 <- 0.25
mu <- 0.25
vals <- rlnorm(100000, mu, sigma2)
hist(vals, breaks=100)

# Variance lnorm
(exp(sigma2) - 1)*exp(2*mu + sigma2)
exp(mu + sigma2/2)
mean(log(vals))

vals2 <- rlnorm(100000, log(95), 0.07332)
hist(vals2, breaks=100)
sqrt(var(vals2))
hist(log(vals2), breaks = 100)

# Is this right?
mu=log(91)
sigma2 = 0.00545
sigma2 = 0.0779
sigma2g = 0.2126

# For the lognormal distribution
exp(mu + sigma2/2) #mean
sqrt((exp(sigma2)-1)*exp(2*mu + sigma2)) #variance

exp(mu + sigma2g/2)
sqrt((exp(sigma2g)-1)*exp(2*mu + sigma2g))

sigmas2 <- seq(0.001, 0.3, by=0.0001)
lnormsd <- sqrt((exp(sigmas2)-1)*exp(2*mu + sigmas2))
ggplot(data.frame(lnormsd=lnormsd,
                  sigmas2 = sigmas2), 
       aes(y=lnormsd,x=sigmas2)) + geom_line()

# Poisson random numbers

# pattern <- c(      1, fung, h1, h2,ins, pred)
# stems <- rpois(100000, 30)
# hist(stems, breaks=100)
# gardsd <- rpois(10000, 5)

spsdg <- 5 # standared deviation for stems between blocks
spsd <- 5  #standard deviation of the residual errors

# What are these, what does 
params <-  c(log(30),    0,  0, -0.1,  0,    0)
treat <- c("c","f","h1","h2","i","p")
gard <- c("g1","g2","g3","g4","g5","g6")
dat <- expand.grid(treat = treat,gard=gard)

lambda <- exp(params[1]+params[2]*(dat$treat == "f")+
                params[3]*(dat$treat == "h1") + 
                params[4]*(dat$treat == "h2") + 
                params[5]*(dat$treat ==  "i") + 
                params[6]*(dat$treat ==  "p")) +
  rep(rnorm(6, 0, 5), each = 6)#random effect of block

dat$abu <- rpois(36, lambda)

abuplot <- ggplot(dat, aes(x = treat, y = abu))
abuplot + geom_jitter(width = 0.1)

fitted_model <- glmer(abu~treat+(1|gard), data=dat, family = "poisson")
summary(fitted_model)



########################################################################



genDatNorm <- function(ngard = 6,    # number of blocks/gardens
                       nplot = 6,    
                       mu = log(95),   # global average log[bio]
                       fung = 0,  # effect size for fungicide
                       ins = 0,  # effect size for insecticide
                       h1 = 0,     # effect size herbivory low
                       h2 = -5,     # effect size herbivory high
                       pred = 0,  # effect size predator
                       sdg = 0.07332, # standrad dev for the random effect
                       sd = 1){    # standard dev for the residuals 
  gard = rep(LETTERS[1:ngard], each = nplot)         # garden label
  gardeff = rep(rnorm(ngard, 0, sdg), each = nplot) # garden effect
  ploteff = rnorm(ngard*nplot, 0, sd)                # noise for plots
  # print(paste("Generating data with sd garden:",sdg," and error var", sd))
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

randDat <- function(rands=9, pattern=c(1, 0, 0,-0.05,0,0), mu=log(95), sdg=0.07332, sd=1, ngard=6){
  # This function uses genDatNorm function to create random dataset and performs linear mixed model analysis. It reports the proportion of datasets where a specified pattern of significance has occured: pattern <- c(int,fung,h1,h2,ins,pred)
  fits <- c()
  dss <- data.frame()
  for (rand in 1:rands){
    ds <- genDatNorm(ngard = ngard, 
                     mu=mu,
                     fung = pattern[2],
                     h1 = pattern[3],
                     h2 = pattern[4],
                     ins=pattern[5],
                     pred = pattern[6],
                     sdg = sdg,
                     sd = sd)
    dss <- rbind(dss, as.data.frame(ds))
    fitted_model <- lmer(bio~treat+(1|gard), data = ds)
    sfm <- summary(fitted_model)
    # print(sfm)
    # print(pattern)
    res <- sum((sfm$coefficients[, "Pr(>|t|)"] < 0.05) == (pattern != 0)) == 6
    fits <- c(fits, res)
  }
  return(list(fits = c(fits), data=dss))
}

mu = 95
mcrands <- 999
sigma2 = 0.0054
sdg = sqrt(sigma2)
sd = sqrt(sigma2)
eff_size = -7
leff_size = log((mu + eff_size)/mu)

res <- randDat(rands=mcrands, 
               pattern=c(1,0,0,leff_size,0,0),
               mu=log(mu),
               sdg,
               sd,
               ngard = 9)
rd <- res$fits


plt <- ggplot(res$data, aes(x = treat, y = exp(bio)))
plt + geom_jitter(width = 0.2, size=2, col="grey50", alpha=0.15) +
  theme_bw() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.2, lwd=1.5) +
  stat_summary(fun.y=mean, geom="point", cex = 5) +
  ylab("Biomass [kg]")

p_hat <- sum(rd/mcrands)
# p_hat + c(-1,1)*1.96*sqrt(p_hat*(1-p_hat)/mcrands)

# Bootstraped confidence intervals
bootSE <- sqrt((sum((rd - mean(rd))^2)/(mcrands -1)))
bootVar <- bootSE^2
p_hat + c(-1,1)*1.96*sqrt(bootVar/mcrands)

#####################################################################

#Curve
mcrands = 999
curveRes <- data.frame()
for(eff_size in -5:-15){
  leff_size = log((mu + eff_size)/mu)
  
  res <- randDat(rands=mcrands, 
                 pattern=c(1,0,0,leff_size,0,0),
                 mu=log(mu),
                 sdg,
                 sd,
                 ngard = 9)
  rd <- res$fits
  # plt <- ggplot(res$data, aes(x = treat, y = exp(bio)))
  # plt + geom_jitter(width = 0.2, size=2, col="grey50", alpha=0.15) +
  #   theme_bw() +
  #   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
  #                geom="errorbar", width=0.2, lwd=1.5) +
  #   stat_summary(fun.y=mean, geom="point", cex = 5) +
  #   ylab("Biomass [kg]")
  
  p_hat <- sum(rd/mcrands)
  # p_hat + c(-1,1)*1.96*sqrt(p_hat*(1-p_hat)/mcrands)
  
  # Bootstraped confidence intervals
  bootSE <- sqrt((sum((rd - mean(rd))^2)/(mcrands -1)))
  bootVar <- bootSE^2
  CIs <- p_hat + c(-1,1)*1.96*sqrt(bootVar/mcrands)
  curveTemp <- data.frame(phat = p_hat, lci = CIs[1], hci = CIs[2])
  curveRes <- rbind(curveRes, curveTemp)
}

# write.table(curveRes, "pCurve_gard9_rands999.txt")

plt <- ggplot(curveRes, aes(x = -15:-5, y = rev(phat)))
plt + geom_line() +
  geom_point() +
  geom_errorbar(aes(ymin=rev(lci), ymax=rev(hci)), width=.2,
              position=position_dodge(.9)) + 
  geom_hline(yintercept=0.8, lty=2, lwd=0.5)




# levelplots ------
library(ggplot2)
library(lmerTest)
library(simr)
library(lattice)
library(latticeExtra)

genDatNorm <- function(ngard = 6,    # number of blocks/gardens
                       nplot = 6,    
                       mu = log(95),   # global average log[bio]
                       fung = 0,  # effect size for fungicide
                       ins = 0,  # effect size for insecticide
                       h1 = 0,     # effect size herbivory low
                       h2 = -5,     # effect size herbivory high
                       pred = 0,  # effect size predator
                       sdg = 0.07332, # standrad dev for the random effect
                       sd = 1){    # standard dev for the residuals 
  gard = rep(LETTERS[1:ngard], each = nplot)         # garden label
  gardeff = rep(rnorm(ngard, 0, sdg), each = nplot) # garden effect
  ploteff = rnorm(ngard*nplot, 0, sd)                # noise for plots
  # print(paste("Generating data with sd garden:",sdg," and error var", sd))
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

randDat <- function(rands=9, pattern=c(1, 0, 0,-0.05,0,0), mu=log(95), sdg=0.07332, sd=1, ngard = 6, nplot = 6){
  # This function uses genDatNorm function to create random dataset and performs linear mixed model analysis. It reports the proportion of datasets where a specified pattern of significance has occured: pattern <- c(int,fung,h1,h2,ins,pred)
  fits <- c()
  dss <- data.frame()
  for (rand in 1:rands){
    ds <- genDatNorm(ngard = ngard,
                     nplot = nplot,
                     mu=mu,
                     fung = pattern[2],
                     h1 = pattern[3],
                     h2 = pattern[4],
                     ins=pattern[5],
                     pred = pattern[6],
                     sdg = sdg,
                     sd = sd)
    dss <- rbind(dss, as.data.frame(ds))
    fitted_model <- lmer(bio~treat+(1|gard), data = ds)
    sfm <- summary(fitted_model)
    # print(sfm)
    # print(pattern)
    res <- sum((sfm$coefficients[, "Pr(>|t|)"] < 0.05) == (pattern != 0)) == 6
    fits <- c(fits, res)
  }
  return(list(fits = c(fits), data=dss))
}

mu=log(95)
sigma2 = 0.0054
sd = sqrt(sigma2)
mu = 95
mcrands <- 99
sdg = sqrt(sigma2)
sd = sqrt(sigma2)
eff_size = -10
leff_size = log((mu + eff_size)/mu)

mcrands <- 99
wbsd <- round(sqrt(seq(0.0005, 0.025, by = 0.005)),3)
bbsd <- round(sqrt(seq(0.0005, 0.025, by = 0.005)),3)
ssize <- seq(5, 20, by=1)#sample size
# leff_size varry the effect size? 
# Random effect standard deviation
resmatgard1 <- matrix(0, 
                      nrow = length(wbsd),
                      ncol = length(ssize))
rownames(resmatgard1) <- as.character(round(wbsd, 3))
colnames(resmatgard1) <- ssize

# ---------------------------------------------------------
temp_res <- randDat(rands=mcrands,
                    pattern=c(1,0,0,leff_size,0,0),
                    mu=log(mu),
                    sdg = bbsd[100],
                    sd = sd,
                    ngard = col,
                    nplot = 6)

plt <- ggplot(temp_res$data, aes(x = treat, y = exp(bio)))
plt + geom_jitter(width = 0.2, size=2, alpha=0.15) +
  theme_bw() +
  ylab("Biomass [kg]")

#----------------------------------------------------------
temp_res <- randDat(rands=mcrands,
                    pattern=c(1,0,0,leff_size,0,0),
                    mu=log(mu),
                    sdg = sdg,
                    sd = wbsd[100],
                    ngard = col,
                    nplot = 6)

plt <- ggplot(temp_res$data, aes(x = treat, y = exp(bio), colour=gard))
plt + geom_jitter(width = 0.2, size=2, alpha=0.15) +
  theme_bw() +
  ylab("Biomass [kg]")



# Sims ----------------------------------------------------
for(row in wbsd){
  for(col in ssize){
    temp_res <- randDat(rands=mcrands,
                        pattern=c(1,0,0,leff_size,0,0),
                        mu=log(mu),
                        sdg = sdg,
                        sd = row,
                        ngard = col,
                        nplot = 6)
    rd <- temp_res$fits
    p_hat <- sum(rd/mcrands)
    
    resmatgard1[as.character(row),as.character(col)] <- p_hat
    
  }
}

# variations ------------------------------------------------
mcrands <- 9
wbsd <- round(sqrt(seq(0.0005, 0.05, by = 0.005)),3)
bbsd <- round(sqrt(seq(0.0005, 0.05, by = 0.005)),3)

resmatgard1 <- matrix(0, 
                      nrow = length(wbsd),
                      ncol = length(bbsd))
rownames(resmatgard1) <- as.character(round(wbsd, 3))
colnames(resmatgard1) <- as.character(round(bbsd, 3))

for(row in wbsd){
  for(col in bbsd){
    temp_res <- randDat(rands=mcrands,
                        pattern=c(1,0,0,leff_size,0,0),
                        mu=log(mu),
                        sdg = col,
                        sd = row,
                        ngard = 6,
                        nplot = 6)
    rd <- temp_res$fits
    p_hat <- sum(rd/mcrands)
    
    resmatgard1[as.character(row),as.character(col)] <- p_hat
    
  }
}

levelplot(resmatgard1, ylab = "Between block variation",
          xlab = "Within block variation")

# write.table(resmatgard1, "bbvsnsize.txt")

levelplot(resmatgard1, xlab = "Between block variation",
          ylab = "Number of blocks")# +
  # as.layer(contourplot(resmatgard1, col = "red"))

for(row in wbsd){
  for(col in ssize){
    temp_res <- randDat(rands=mcrands,
                        pattern=c(1,0,0,leff_size,0,0),
                        mu=log(mu),
                        sdg = sdg,
                        sd = row,
                        ngard = col,
                        nplot = 6)
    rd <- temp_res$fits
    p_hat <- sum(rd/mcrands)
    
    resmatgard1[as.character(row),as.character(col)] <- p_hat
    
  }
}

levelplot(resmatgard1, xlab = "Within block variation",
          ylab = "Number of blocks") +
  as.layer(contourplot(resmatgard1, col = "red"))


plt <- ggplot(temp_res$data, aes(x = treat, y = exp(bio)))
plt + geom_jitter(width = 0.2, size=2, col="grey50", alpha=0.15) +
  theme_bw() +
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
               geom="errorbar", width=0.2, lwd=1.5) +
  stat_summary(fun.y=mean, geom="point", cex = 5) +
  ylab("Biomass [kg]")


res <- genDatNorm(ngard = 20,
           nplot = 6,
           mu = log(mu),   
           fung = 0,  
           ins = 0, 
           h1 = 0,     
           h2 = 5*leff_size,     
           pred = 0,
           sdg = 0.0005,
           sd = 0.5)

res2 <- genDatNorm(ngard = 20,
                   nplot = 6,
                   mu = log(mu),   
                   fung = 0,  
                   ins = 0, 
                   h1 = 0,     
                   h2 = 5*leff_size,     
                   pred = 0,
                   sdg = 0.5,
                   sd = 0.0005)

res$ss <- "black"
res2$ss <- "yellow"

dat <- rbind(res, res2)

summary(lmer(bio~treat+(1|gard), data = res))



plt <- ggplot(dat, aes(x = treat, y = exp(bio), colour = ss))
plt + geom_jitter(width = 0.2, size=2, alpha=0.15) +
  theme_bw() +
  ylab("Biomass [kg]")





# Once again with simr ------------
# Run the analyses_wanang.R script first

library("simr")
# powerSim(bio_rbl)
ds <- getData(bio_rbl)
dsnull <- ds[,c("GARDEN", "TREAT")]

#########################################################################

eff_size = -10
mu = 95
leff_size = log((mu + eff_size)/mu)
treat <- c("c","f","h1","h2","i","p")
gard <- c("g1","g2","g3","g4","g5","g6")


dat <- expand.grid(gard=gard, treat = treat)

expected <- makeLmer(y~treat+(1|gard),
         fixef = c(log(95),0,0,leff_size,0,-leff_size), 
         VarCorr = 0.005,
         sigma = sqrt(0.005),
         data=dat)

#########################################################################
# Change the sample size and the variance

# predicted variance
sigma2 = 0.0054
sdg <- sqrt(sigma2)

# doubling the number of samples
ssize <- 9

# Change fixed effects
eff_size = -15
mu = 95

leff_size = log((mu + eff_size)/mu)

fixef(bio_rbl) <- c(log(95),0,0,0,leff_size,0)

# VarCorr(bio_rbl)# <- matrix(c(1,0.5,0.5,1), nrow=2, byrow = T)
# VarCorr(bio_rbl) <- matrix(0.5, nrow=2, ncol=2) + diag(0.5, 2)

VarCorr(bio_rbl) <- 0.0052
sigma(bio_rbl) <- sqrt(0.0052)

# Create the dataset
vals <- doSim(bio_rbl)
dat <- cbind(vals,  dsnull)
plt <- ggplot(dat, aes(x = TREAT, y = exp(vals), colour = GARDEN))
plt + geom_jitter(width = 0.2, size=2, alpha=0.15) +
  theme_bw() +
  ylab("Biomass [kg]")

powerSim(bio_rbl)

vars <- seq(0.003, 0.090, by=0.015)
lognormSD <- sqrt(exp(sqrt(vars) - 1) * exp(2*mu+sqrt(vars)))


resmatgard1 <- matrix(0, 
                      nrow = length(vars),
                      ncol = length(vars))
rownames(resmatgard1) <- as.character(round(vars, 3))
colnames(resmatgard1) <- as.character(round(vars, 3))

mcrands <- 99
counter <- 1
dims    <- length(vars)
for (i in vars){
  for(j in vars){
    counter <- counter + 1
    print(counter/(dims*dims))
    VarCorr(bio_rbl) <- i  #
    sigma(bio_rbl) <- sqrt(j)
    est <- powerSim(bio_rbl, nsim = mcrands)
    resmatgard1[as.character(i),
                as.character(j)] <- summary(est)$mean
  }
}

levelplot(resmatgard1, ylab = "Sigma",
          xlab = "VarCorr")


bio_ext_gardens <- extend(bio_rbl, within="GARDEN", n=20)
pC <- powerCurve(bio_ext_gardens)
plot(pC)




##### STOLEN FROM THE ORIGINAL
# 
# # Loading packages and defining functions
# 
# library(ggplot2)
# library(lmerTest)
# library(simr)
# library(lattice)
# library(latticeExtra)
# 
# genDatNorm <- function(ngard = 6,    # number of blocks/gardens
#                     nplot = 6,    
#                     mu = log(95),   # global average log[bio]
#                     fung = 0,  # effect size for fungicide
#                     ins = 0,  # effect size for insecticide
#                     h1 = 0,     # effect size herbivory low
#                     h2 = -5,     # effect size herbivory high
#                     pred = 0,  # effect size predator
#                     sdg = 0.07332, # standrad dev for the random effect
#                     sd = 1){    # standard dev for the residuals 
#   gard = rep(LETTERS[1:ngard], each = nplot)         # garden label
#   gardeff = rep(rnorm(ngard, 0, sdg), each = nplot) # garden effect
#   ploteff = rnorm(ngard*nplot, 0, sd)                # noise for plots
#   # print(paste("Generating data with sd garden:",sdg," and error var", sd))
#   # Simulate the explanatory variables
#   treat <- c()
#   for (g in 1:ngard){treat <- c(treat,
#                                 sample(c("c","f","i","h1","h2","p"), nplot))}
#   insecticide <- as.numeric(treat == "i")
#   fungicide   <- as.numeric(treat == "f")
#   herbivory1  <- as.numeric(treat == "h1")
#   herbivory2  <- as.numeric(treat == "h2")
#   predator    <- as.numeric(treat == "p")
#   bio <- mu + fung*fungicide + ins*insecticide + h1*herbivory1 + h2*herbivory2 + pred*predator + gardeff + ploteff
#   ds <- as.data.frame(cbind(bio,
#                             gard,
#                             treat))
#   
#   ds$bio <- as.numeric(as.character(ds$bio))
#   return(ds)
# }
# 
# randDat <- function(rands=9, pattern=c(1, 0, 0,-0.05,0,0), mu=log(95), sdg=0.07332, sd=1, ngard = 6, nplot = 6){
#   # This function uses genDatNorm function to create random dataset and performs linear mixed model analysis. It reports the proportion of datasets where a specified pattern of significance has occured: pattern <- c(int,fung,h1,h2,ins,pred)
#   fits <- c()
#   dss <- data.frame()
#   for (rand in 1:rands){
#     ds <- genDatNorm(ngard = ngard,
#                      nplot = nplot,
#                      mu=mu,
#                      fung = pattern[2],
#                      h1 = pattern[3],
#                      h2 = pattern[4],
#                      ins=pattern[5],
#                      pred = pattern[6],
#                      sdg = sdg,
#                      sd = sd)
#     dss <- rbind(dss, as.data.frame(ds))
#     fitted_model <- lmer(bio~treat+(1|gard), data = ds)
#     sfm <- summary(fitted_model)
#     # print(sfm)
#     # print(pattern)
#     res <- sum((sfm$coefficients[, "Pr(>|t|)"] < 0.05) == (pattern != 0)) == 6
#     fits <- c(fits, res)
#   }
#   return(list(fits = c(fits), data=dss))
# }



# res <- randDat(rands=mcrands, 
#                pattern=c(1,0,0,leff_size,0,0),
#                mu=log(mu),
#                sdg,
#                sd)
# rd <- res$fits
# 
# #An example dataset simulated based on our variance estimates.
# # An example simulat
# # ds <- genDat(ins=-10)
# plt <- ggplot(res$data, aes(x = treat, y = exp(bio)))
# plt + geom_jitter(width = 0.2, size=2, col="grey50", alpha=0.15) +
#   theme_bw() +
#   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#                  geom="errorbar", width=0.2, lwd=1.5) +
#   stat_summary(fun.y=mean, geom="point", cex = 5) +
#   ylab("Biomass [kg]")
# 
# p_hat <- sum(rd/mcrands)
# # p_hat + c(-1,1)*1.96*sqrt(p_hat*(1-p_hat)/mcrands)
# 
# # Bootstraped confidence intervals
# bootSE <- sqrt((sum((rd - mean(rd))^2)/(mcrands -1)))
# bootVar <- bootSE^2
# # p_hat + c(-1,1)*1.96*sqrt(bootVar/mcrands)


# # Poisson random numbers
# spsdg <- 5 # standared deviation for stems between blocks
# spsd <- 5  #standard deviation of the residual errors
# 
# # What are these, what does 
# params <-  c(log(30),    0,  0, -0.1,  0,    0)
# 
# treat <- c("c","f","h1","h2","i","p")
# gard <- c("g1","g2","g3","g4","g5","g6")
# dat <- expand.grid(treat = treat,gard=gard)
# 
# datdf <- data.frame()
# for (rand in 1:mcrands){
#   
#   lambda <- exp(params[1]+params[2]*(dat$treat == "f")+
#                 params[3]*(dat$treat == "h1") + 
#                 params[4]*(dat$treat == "h2") + 
#                 params[5]*(dat$treat ==  "i") + 
#                 params[6]*(dat$treat ==  "p")) +
#   rep(rnorm(6, 0, 5), each = 6)#random effect of block
#   dat$abu <- rpois(36, lambda)
#   datdf <- rbind(datdf, dat)
# }
# 
# 
# lambda <- exp(params[1]+params[2]*(dat$treat == "f")+
#                 params[3]*(dat$treat == "h1") + 
#                 params[4]*(dat$treat == "h2") + 
#                 params[5]*(dat$treat ==  "i") + 
#                 params[6]*(dat$treat ==  "p")) +
#   rep(rnorm(6, 0, 5), each = 6)#random effect of block
# dat$abu <- rpois(36, lambda)
# 
# abuplot <- ggplot(datdf, aes(x = treat, y = abu))
# abuplot + geom_jitter(width = 0.2, alpha=0.1, col="grey50") + 
#   xlab("Treatment") + ylab("No of stems per 25 m^2 plot")
# 
# #fitted_model <- glmer(abu~treat+(1|gard), data=dat, family = "poisson")
# #summary(fitted_model)




library(simr)
library(ggplot2)

eff_size = -10
mu = 95
leff_size = log((mu + eff_size)/mu)
treat <- c("c","f","h1","h2","i","p")
gard <- LETTERS[1:20]

dat <- expand.grid(gard=gard, treat = treat)

bioexp <- makeLmer(y~treat+(1|gard),
                   fixef = c(log(95),0,0,leff_size,-leff_size,-leff_size/2), 
                   VarCorr = 0.01,
                   sigma = sqrt(0.01),
                   data=dat)


dt <- getData(bioexp)

plt <- ggplot(dt, aes(x = treat, y = y, colour = gard))
plt + geom_jitter(width = 0.2, size=2, alpha=0.15) +
  theme_bw() +
  ylab("Biomass [kg]")

pc2 <- powerCurve(bioexp, along = "gard", breaks=c(3,6,9,12), nsim=99)
plot(pc2)


mu = 30
mcrands <- 99
sdg = 5
sd = 5
eff_size = -5
leff_size = log((mu + eff_size)/mu)

# Create dataset
treat <- c("c","f","h1","h2","i","p")
gard <- LETTERS[1:20]

dat <- expand.grid(gard=gard, treat = treat)

res <- c()
for(i in 1:100){
  richexp <- makeGlmer(y~treat+(1|gard),
                       family = "poisson",
                       fixef = c(log(mu),0,0,
                                 leff_size,-leff_size,-leff_size/2), 
                       VarCorr = 0.0225,
                       data=dat)
  
  dt <- getData(richexp)
  res <- c(res, sd(tapply(dt$y, dt$gard, mean)))
}

mean(res)

plt <- ggplot(dt, aes(x = treat, y = y, colour = gard))
plt + geom_jitter(width = 0.2, size=2, alpha=0.15) +
  theme_bw() +
  ylab("Species richness")





tagb <- 100
mu=log(tagb)
sigma2 = 0.01#20
sd = sqrt(sigma2)

# Histogram of biomasses
vals2 <- rlnorm(10000, mu, sd) # inputs are the mu and sd of the noorm vals
hist(vals2, breaks=100, probability=T)
sd(vals2)


# Levelplot dla rozkladu poissona

mu = 30
sdg = 6
sd = 6
eff_size = 6
leff_size = log((mu + eff_size)/mu)

# Create dataset
treat <- c("c","f","h1","h2","i","p")
gard <- LETTERS[1:20]

dat <- expand.grid(gard=gard, treat = treat)

richexp <- makeGlmer(y~treat+(1|gard),
                     family = "poisson",
                     fixef = c(log(mu),                    # control
                               -leff_size,                # fungi
                               (leff_size/2+leff_size)/2, # herb moderate
                               leff_size,                 # herbivory high
                               -leff_size,                # insecticide
                               leff_size/2),              # predator, 
                     VarCorr = 0.0225, # estimated such that s_e is ~ 5. 0.0225
                     data=dat)

dt <- getData(richexp)

par(mfrow=c(1,2))
plt <- ggplot(dt, aes(x = treat, y = y, colour = gard))
plt + geom_jitter(width = 0.2, size=2, alpha=0.5) +
  theme_bw() +
  ylab("Species richness") + 
  xlab("Treatment")

mus <- seq(2,10,by=1)

respois <- c()
for(mu in mus){
  print(mu)
  fixef(richexp) = c(log(mu),                    # control
            -leff_size,                # fungi
            (leff_size/2+leff_size)/2, # herb moderate
            leff_size,                 # herbivory high
            -leff_size,                # insecticide
            leff_size/2)
  est <- powerSim(richexp, nsim = mcrands, progress = F)
  respois <- c(respois,summary(est)$mean)
}

respois

# Rename rows and columns
newnames <- c()
for(value in vars){
  vals <- rlnorm(10000, log(95), sqrt(value))
  newnames <- c(newnames, sd(vals))
}

newnames <- round(newnames, 0)
rownames(resmatgard1) <- as.character(newnames)
colnames(resmatgard1) <- as.character(newnames)

levelplot(resmatgard1, ylab = "Within block variation (Sigma)",
          xlab = "Between block variation (VarCorr)") +
  as.layer(contourplot(resmatgard1, col = "red"))
