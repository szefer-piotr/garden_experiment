# Specify a path for data and rds files
setwd("C:\\Users\\Piotr Szefer\\Desktop\\Work\\garden experiment\\figs\\stanplots")

# Examples
# https://www.sciencedirect.com/science/article/pii/S0095447017302310

# Upload data
AllTestData <- read.table("all_test_data.txt")
TreeTestDataset <- read.table("tree_test_data.txt")

# Load package 
library(brms)

# Perform the tests only one example for log(Bio)

## Set some priors 
priors <- c(set_prior("normal(0, 200)", class = "Intercept"),
            set_prior("normal(0, 50)", class = "b"),
            set_prior("normal(0, 100)", class = "sd"),
            set_prior("normal(0, 100)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

## Whole community random intercept
l1W <- brm(log(BIO)~TREAT + (1|GARDEN), 
           data=AllTestData,
           #prior = priors,
           control = list(adapt_delta = 0.999),
           file="bioAll")

## Whole community random intercept and slope
l1W_ranef <- brm(log(BIO) ~ TREAT + (1+TREAT|GARDEN), 
           data=AllTestData,
           #prior = priors,
           control = list(adapt_delta = 0.999),
           file="bioAllrenef")

l1W_ranslope <- brm(log(BIO) ~ TREAT + (0+TREAT|GARDEN), 
                 data=AllTestData,
                 #prior = priors,
                 control = list(adapt_delta = 0.999),
                 file="bioAllrenslope")

# Summaries
summary(l1W)
summary(l1W_ranef)

# See if the model with random intercepts is more informative than simple mixed model
compare_ic(waic(l1W_ranef), waic(l1W_ranslope), ic = "waic") # it looks like it is!

# Plots
stanplot(l1W, type="hist")                    # random intercept
stanplot(l1W_ranef, pars="^b_", type = "hist")  # random intercept and slope
marginal_effects(l1W)      
marginal_effects(l1W_ranef)      # second model seems to better represent the data

# Creating plots for the predictions

# new data for the model with random intercept
newdata1 = data.frame(TREAT = levels(AllTestData$TREAT))

fit1 = fitted(
  l1W,
  newdata=newdata1,
  re_formula = NA, # ignore random effects
  summary = TRUE   # mean and 95% ci
)

colnames(fit1) = c('fit', 'se','lwr','upr')
df_plot1 = cbind(newdata1, fit1)


# new data for the model with random intercept and slope
newdata = expand.grid(TREAT= levels(AllTestData$TREAT),
                      GARDEN = levels(AllTestData$GARDEN))

fit = fitted(
  l1W_ranef,
  newdata=newdata1,
  re_formula = NA, # ignore random effects
  summary = TRUE   # mean and 95% ci
)

# fit = fitted(
#   l1W,
#   newdata=newdata,
#   re_formula = NULL, # ignore random effects
#   summary = TRUE   # mean and 95% ci
# )

colnames(fit) = c('fit', 'se','lwr','upr')
df_plot = cbind(newdata, fit)

# Averaging different garden means and CIs for plotting
df_mean = tapply(df_plot$fit, df_plot$TREAT, mean)
df_lwr = tapply(df_plot$lwr, df_plot$TREAT, mean)
df_upr = tapply(df_plot$upr, df_plot$TREAT, mean)
df_collapsed = cbind(df_mean, df_lwr, df_upr)
colnames(df_collapsed) = c('fit','lwr','upr')
df_collapsed <- as.data.frame(df_collapsed)
df_collapsed$TREAT = rownames(df_collapsed)

ggplot(df_collapsed, aes(x = TREAT, y = fit)) +
  #geom_violin(data=AllTestData, aes(x=TREAT, y=log(BIO)), alpha=0.5, color="gray70", fill='gray95') +
  geom_jitter(data=AllTestData, aes(x=TREAT, y=log(BIO)), alpha=0.3, size=4,
              position = position_jitter(width = 0.07)) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position=position_dodge(), size=1, width=.5) +
  geom_point(shape=21, size=4, fill='red') +
  xlab("") +
  ylab('log(BIO)') +
  theme_bw () +
  theme(panel.grid = element_blank())

#
# ggplot(df_plot, aes(x = TREAT, y = fit)) +
#   #geom_violin(data=AllTestData, aes(x=TREAT, y=log(BIO)), alpha=0.5, color="gray70", fill='gray95') +
#   geom_jitter(data=AllTestData, aes(x=TREAT, y=log(BIO)), alpha=0.3, size=4,
#               position = position_jitter(width = 0.07)) +
#   geom_errorbar(aes(ymin=lwr, ymax=upr), position=position_dodge(), size=1, width=.5) +
#   geom_point(shape=21, size=4, fill='red') +
#   xlab("") +
#   ylab('log(BIO)') +
#   theme_bw () +
#   theme(panel.grid = element_blank())


###############################
# Facet plots for each garden #
###############################


# First for the model with random intercept

# new data for the model with random intercept
newdata = expand.grid(TREAT= levels(AllTestData$TREAT),
                      GARDEN = levels(AllTestData$GARDEN))

# fit = fitted(
#   l1W_ranef,
#   newdata=newdata,
#   re_formula = NULL, # ignore random effects
#   summary = TRUE   # mean and 95% ci
# )

fit = fitted(
  l1W,
  newdata=newdata,
  re_formula = NULL, # ignore random effects
  summary = TRUE   # mean and 95% ci
)

colnames(fit) = c('fit', 'se','lwr','upr')
df_plot = cbind(newdata, fit)

ggplot(df_plot, aes(x = TREAT, y = fit)) +
  #geom_violin(data=AllTestData, aes(x=TREAT, y=log(BIO)), alpha=0.5, color="gray70", fill='gray95') +
  geom_jitter(data=AllTestData, aes(x=TREAT, y=log(BIO)), alpha=0.3, size=4,
              position = position_jitter(width = 0.07)) +
  geom_point(shape=21, size=4, fill='red') +
  facet_wrap(~GARDEN, scales="free", drop=TRUE)+
  xlab("") +
  ylab('log(BIO)') +
  theme_bw () +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle=90,hjust=1, size=7))


# And second for the model with random intercept

# new data for the model with random intercept
newdata = expand.grid(TREAT= levels(AllTestData$TREAT),
                      GARDEN = levels(AllTestData$GARDEN))

fit = fitted(
  l1W_ranef,
  newdata=newdata,
  re_formula = NULL, # ignore random effects
  summary = TRUE   # mean and 95% ci
)

# fit = fitted(
#   l1W,
#   newdata=newdata,
#   re_formula = NULL, # ignore random effects
#   summary = TRUE   # mean and 95% ci
# )

colnames(fit) = c('fit', 'se','lwr','upr')
df_plot = cbind(newdata, fit)

ggplot(df_plot, aes(x = TREAT, y = fit)) +
  #geom_violin(data=AllTestData, aes(x=TREAT, y=log(BIO)), alpha=0.5, color="gray70", fill='gray95') +
  geom_jitter(data=AllTestData, aes(x=TREAT, y=log(BIO)), alpha=0.3, size=4,
              position = position_jitter(width = 0.07)) +
  geom_point(shape=21, size=4, fill='red') +
  facet_wrap(~GARDEN, scales="free", drop=TRUE)+
  xlab("") +
  ylab('log(BIO)') +
  theme_bw () +
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle=90,hjust=1, size=7))


#############
# Contrasts #
#############

# ONE EXAMPLE (out of 45) of a pairwise comparison.

# Sampling differences, posterior for the differences between C and I
fit1 = as.data.frame(fitted(l1W_ranef,
                            newdata=newdata,
                            re_formula = NULL, # include all effects?
                            summary = FALSE #extract the full MCMC
))
colnames(fit1) <- newdata$TREAT

# Regardless of the random effects of the block and not allowing for different responces to 
# treatments in different blocks.
#fit1 = as.data.frame(fitted(l1W_ranef,
#                            newdata=newdata,
#                            re_formula = NULL, # include all effects?
#                            summary = FALSE #extract the full MCMC
#))
#colnames(fit1) <- newdata$TREAT

# Is CONTROL different from PREDATOR regardless of the block?
covspre <- fit1[ ,colnames(fit1) == "CONTROL"] - fit1[ ,colnames(fit1) == "PREDATOR"]

# Histogram for the differences
hist(apply(covspre, 1, mean), breaks = 50,
     main = "Posterior: CONTROL and PREDATOR tratments",
     xlab="Difference")

# Probability of observing difference between Control and Predator above zero
# from different gardens
sum(covspre>0)/(6*4000) # 86.6 % probablity that it is so


## NOTES ##
    ###
    ###
  #######
   #####
    ###
     #

l2W <- brm(SW ~ TREAT,
                 data=AllTestData,
                 control = list(adapt_delta = 0.999),
                 file="shaAll2")

l2W_ranef <- brm(SW ~ TREAT + (1+TREAT|GARDEN),
           data=AllTestData,
           control = list(adapt_delta = 0.97),
           file="shaAll_ranef")

compare_ic(waic(l2W), waic(l2W_ranef), ic = "waic")

summary(l2W)
summary(l2W_ranef)
stanplot(l2W, pars="^b_")

swnewdata <- data.frame(TREAT= levels(AllTestData$TREAT))

swfit = fitted(
  l2W_ranef,
  newdata=swnewdata,
  re_formula = NA, # ignore random effects
  summary = TRUE   # mean and 95% ci
)

colnames(swfit) = c('fit', 'se','lwr','upr')
swdf_plot = cbind(swnewdata, swfit)

ggplot(swdf_plot, aes(x = TREAT, y = fit)) +
  #geom_violin(data=AllTestData, aes(x=TREAT, y=SW), alpha=0.5, color="gray70", fill='gray95') +
  geom_jitter(data=AllTestData, aes(x=TREAT, y=SW), size=4,
              alpha=0.3, position = position_jitter(width = 0.07)) +
  geom_errorbar(aes(ymin=lwr, ymax=upr), position=position_dodge(), size=1, width=.5) +
  geom_point(shape=21, size=4, fill='red') +
  xlab("") +
  ylab('SW') +
  theme_bw () +
  theme(panel.grid = element_blank())

 
# l3W <- brm(SPEC_NO~TREAT+(1|GARDEN), data=AllTestData,
#            family = poisson(),
#            control = list(adapt_delta = 0.99),
#            file="nosAll")
# 
# l4W <- brm(EVEN~TREAT+(1|GARDEN), data=AllTestData,
#            family = Beta(),
#            control = list(adapt_delta = 0.99),
#            file="eveAll")
# 
# #### Only trees
# l1Wt <- brm(logBio~TREATMENT+(1|GARDEN), data=TreeTestDataset,
#            control = list(adapt_delta = 0.99),
#            file="bioTre")
# 
# l2Wt <- brm(DivTree~TREATMENT + (1|GARDEN), data=TreeTestDataset,
#             control = list(adapt_delta = 0.99),
#             file="swiTre")
# 
# l3Wt <- brm(species~TREATMENT+(1|GARDEN), data=TreeTestDataset,
#               family = poisson(),
#               control = list(adapt_delta = 0.99),
#               file="nosTre")
# 
# l4Wt <- brm(Evenness~TREATMENT + (1|GARDEN), data=TreeTestDataset,
#              family = Beta(),
#              control = list(adapt_delta = 0.99),
#              file = "eveTre")
# 
# l5Wt <- brm(stems~TREATMENT + (1|GARDEN), data=TreeTestDataset,
#             family=poisson(),
#             control = list(adapt_delta = 0.9999),
#             file="steTre")
# 
# 
# bioAll <- stanplot(l1W, pars="^b_")
# swiAll <- stanplot(l2W, pars="^b_")
# nosAll <- stanplot(l3W, pars="^b_")
# eveAll <- stanplot(l4W, pars="^b_")
# bioTre <- stanplot(l1Wt, pars="^b_")
# swiTre <- stanplot(l2Wt, pars="^b_")
# nosTre <- stanplot(l3Wt, pars="^b_")
# eveTre = stanplot(l4Wt, pars="^b_")
# steTre = stanplot(l5Wt, pars="^b_")