

# Examples
# https://www.sciencedirect.com/science/article/pii/S0095447017302310


# load data and packages --------------------------------------------------

# load data
AllTestData <- read.table("figs/stanplots/all_test_data.txt")
TreeTestDataset <- read.table("figs/stanplots/tree_test_data.txt")

# Load package 
library(brms)
library(ggplot2)
library(bayesplot)
library(tidybayes)
library(dplyr)


# plots of the raw data ---------------------------------------------------

ggplot(AllTestData, aes(y = TREAT, x = BIO)) + 
  geom_point() + 
  facet_wrap(~GARDEN)


# setting various priors --------------------------------------------------

priors_raninter <- c(set_prior("normal(0, 200)", class = "Intercept"),
                     set_prior("normal(0, 50)", class = "b"),
                     set_prior("normal(0, 5)", class = "sd"),
                     set_prior("normal(0, 5)", class = "sigma"))

priors_rsi <- c(set_prior("normal(0, 200)", class = "Intercept"),
            set_prior("normal(0, 50)", class = "b"),
            set_prior("normal(0, 5)", class = "sd"),
            set_prior("normal(0, 5)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

priors_raninter2 <- c(set_prior("normal(4, 2)", class = "Intercept"),
                 set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 4)", class = "sd"),
                 set_prior("normal(0, 3)", class = "sigma"))

priors_rsi2 <- c(set_prior("normal(4, 2)", class = "Intercept"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 4)", class = "sd"),
            set_prior("normal(0, 3)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))


# ANDREW random intercepts and slopes --------------------------------------------

#define model 
#does changing contrasts change the model?
levels(AllTestData$TREAT)

#contrasts(AllTestData$TREAT) <- contr.helmert(n = 6)

bio_lnorm_bf <- bf(BIO ~ 1 + TREAT + (1 + TREAT|GARDEN), family = lognormal())
brms::get_prior(bio_lnorm_bf, data = AllTestData)

## Set some priors 
priors <- c(set_prior("normal(4, 2)", class = "Intercept"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 4)", class = "sd"),
            set_prior("normal(0, 3)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

## Whole community with random intercept and slope
l1W_ranef_lnorm <- brm(bio_lnorm_bf, 
                 data=AllTestData,
                 prior = priors_rsi2,
                 control = list(adapt_delta = 0.999),
                 file="bioAllrenef_lnormal", cores = 4, sample_prior = TRUE)

stancode(l1W_ranef)

#what are the parameters here worth thinking about?
brms::parnames(l1W_ranef)
#get the info on divergent iterations
np <- nuts_params(l1W_ranef)
# customize appearance of divergences
color_scheme_set("darkgray")
div_style <- parcoord_style_np(div_color = "orange", div_size = 1, div_alpha = 1)
mcmc_parcoord(as.array(l1W_ranef),
              size = 0.25, alpha = 0.1,
              regex_pars = "sd_.*",
              np = np, np_style = div_style)+ coord_flip()

l1W_ranef %>% summary

# recoding factors --------------------------------------------------------

#the model above still has a hard time with divergent iterations. I am curious
#to try how it works if the variables are coded by hand to develop contrasts

AllTestData_recode <- AllTestData %>% 
  #remove the weevil stuff just for now
  filter(!(TREAT %>% str_detect("WEEVIL"))) %>% 
  # removal treatments as several columns
  mutate(enemy_removal = if_else(TREAT %in% c("FUNGICIDE", "INSECTICIDE", "PREDATOR"), true = 1, false = 0),
         is_spray      = if_else(TREAT %in% c("FUNGICIDE", "INSECTICIDE"),             true = 1, false = 0),
         insecticide   = if_else(TREAT %in% c("INSECTICIDE"),                          true = 1, false = 0))



bio_lnorm_recode_bf <- bf(BIO ~ 1 + enemy_removal + is_spray + insecticide + (1 + enemy_removal + is_spray + insecticide |GARDEN), family = lognormal())
brms::get_prior(bio_lnorm_recode_bf, data = AllTestData_recode)
## Set some priors 
priors <- c(set_prior("normal(4, 2)", class = "Intercept"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 2)", class = "sd"),
            set_prior("normal(0, 3)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

## Whole community with random intercept and slope
bio_recode <- brm(bio_lnorm_recode_bf, 
                 data=AllTestData_recode,
                 prior = priors,
                 control = list(adapt_delta = 0.9),
                 file="bio_recode", cores = 4, sample_prior = TRUE)

mcmc_parcoord(as.array(bio_recode),
              size = 0.25, alpha = 0.1,
              regex_pars = "sd_.*",
              np = nuts_params(bio_recode), np_style = div_style)+ coord_flip()

# experiment with gamma distribution for biomass instead ------------------


#perhaps not necessary, since using lognormal distribution seems to make a big difference.


# weevils as a number ------------------------------------------------------

AllTestData_numweevil <- AllTestData %>% 
  # could use relevel here
  mutate(removal_trt = if_else(TREAT %in% c("CONTROL", "WEEVIL125", "WEEVIL25"), true = "0_removal", false = as.character(TREAT)),
         removal_trt = as.factor(removal_trt),
         weevil = readr::parse_number(TREAT)) %>% 
  tidyr::replace_na(list(weevil = 0)) %>% 
  mutate(weevil = weevil / 25)

bio_lnorm_numweev_bf <- bf(BIO ~ 1 + removal_trt + weevil + (1 + removal_trt + weevil|GARDEN), family = lognormal())
brms::get_prior(bio_lnorm_numweev_bf, data = AllTestData_numweevil)
## Set some priors 
priors <- c(set_prior("normal(4, 2)", class = "Intercept"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 4)", class = "sd"),
            set_prior("normal(0, 3)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

## Whole community with random intercept and slope
l1W_ranef_numweev <- brm(bio_lnorm_numweev_bf, 
                 data=AllTestData_numweevil,
                 prior = priors,
                 control = list(adapt_delta = 0.95),
                 file="l1W_ranef_numweev", cores = 4, sample_prior = TRUE)

summary(l1W_ranef_numweev)

marginal_effects(l1W_ranef_numweev)

# ic comparison -----------------------------------------------------------



# See if the model with random intercepts is more informative than simple mixed model
compare_ic(waic(l1W), waic(l1W_ranef), ic = "waic") # it looks like it is!



# PIOTR -----------------------------------------------------------------------


## Biomass -----------------------------------------------------

# !!! Using more informative priors

# Random intercept
l1W <- brm(log(BIO)~TREAT + (1|GARDEN), 
           data=AllTestData,
           prior = priors_raninter2,
           control = list(adapt_delta = 0.999),
           file="bioAll")

## Whole community with random intercept and slope
l1W_ranef <- brm(log(BIO) ~ TREAT + (1+TREAT|GARDEN), 
                 data=AllTestData,
                 prior = priors_rsi2,
                 control = list(adapt_delta = 0.999),
                 file="bioAllrenef")

# Comparisons

compare_ic(waic(l1W), waic(l1W_ranef), ic = "waic")

# model diagnostics -------------------------------------------------------------

pp = brms::pp_check(l1W_ranef)
pp + theme_bw()


# plots -------------------------------------------------------------------

ranef_plot <- marginal_effects(l1W_ranef)    # second model seems to better represent the data

# Creating plots for the predictions -------------------------------------

# Using predictions from the stanplots they produce the same results as predictions
ranef_plot$TREAT

## for comparison another plot using predictions of marginal_effect plot

# Plot template
ggplot(ranef_plot$TREAT, aes(x = TREAT, y = estimate__)) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), position=position_dodge(), size=1, width=.5) +
  geom_jitter(data=AllTestData, aes(x=TREAT, y=log(BIO), group=GARDEN), 
              alpha=0.10, size=4,
              position = position_jitter(width = 0.07)) +
  geom_line(data=AllTestData, aes(x=TREAT, y=log(BIO), 
                                  group=GARDEN, 
                                  linetype = GARDEN), 
            size = 0.5, alpha=0.25) +
  geom_point(shape=21, size=4, fill='red') +
  xlab("") +
  ylab('log(BIO)') +
  theme_bw () +
  theme(panel.grid = element_blank())

# Pairwise comparisons -----------------------------------------------------------

sph <- stanplot(l1W_ranef, type = "hist", pars="^b_")

par(mfrow = c(3,2))

# Change the name of this
sph$data$Parameter

for(type in unique(sph$data$Parameter)){
  
  vs <- sph$data$value[sph$data$Parameter == type]
  
  colour = "lightblue"
  
  if (mean(vs<0) > 0.90) colour = "orange"
  if (mean(vs<0) > 0.95) colour = "red"
  
  hist(vs, breaks = 50,
       main = paste(type, " ", round(mean(vs<0)*100,1), " %"),
       xlab="Difference", col = colour)
  
  abline(v = 0, lty=2)
  
}

stanplot(l1W_ranef, type = "hist",  pars="^b_")

# Histogram for the differences

# Rest of the models --------------------------------------------------------

# SHANNON ------------------------------------------------------------------

sw_norm_f <- bf(SW ~ 1 + TREAT + (1 + TREAT|GARDEN), family = gaussian())
sw_prior <- get_prior(sw_norm_f, data = AllTestData)
sw_prior$prior[1] <- "normal(1.25,4)"

# priors <- c(set_prior("normal(4, 2)", class = "Intercept"),
#             set_prior("normal(0, 3)", class = "b"),
#             set_prior("normal(0, 4)", class = "sd"),
#             set_prior("normal(0, 3)", class = "sigma"),
#             set_prior("lkj(2)", class = "cor"))

# How to get better 
l2W_ranef <- brm(SW ~ TREAT + (1+TREAT|GARDEN),
           data=AllTestData,
           control = list(adapt_delta = 0.999),
           prior = sw_prior,
           file = "swAllranef")


#plot(l2W_ranef)

pp = brms::pp_check(l2W_ranef)
pp + theme_bw()

sw_ranef_plot <- marginal_effects(l2W_ranef)    # second model seems to better represent the data

## for comparison another plot using predictions of marginal_effect plot
# Plot template
ggplot(sw_ranef_plot$TREAT, aes(x = TREAT, y = estimate__)) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), position=position_dodge(), size=1, width=.5) +
  geom_jitter(data=AllTestData, aes(x=TREAT, y=SW, group=GARDEN), 
              alpha=0.10, size=4,
              position = position_jitter(width = 0.07)) +
  geom_line(data=AllTestData, aes(x=TREAT, y=SW, 
                                  group=GARDEN, 
                                  linetype = GARDEN), 
            size = 0.5, alpha=0.25) +
  geom_point(shape=21, size=4, fill='red') +
  xlab("") +
  ylab("Shannon's Index") +
  theme_bw () +
  theme(panel.grid = element_blank())

# Pairwise comparisons

sph <- stanplot(l2W_ranef, type = "hist", pars="^b_")

par(mfrow = c(3,2))

# Change the name of this
sph$data$Parameter

for(type in unique(sph$data$Parameter)){
  
  vs <- sph$data$value[sph$data$Parameter == type]
  
  colour = "lightblue"
  
  if (mean(vs<0) > 0.90) colour = "orange"
  if (mean(vs<0) > 0.95) colour = "red"
  
  hist(vs, breaks = 50,
       main = paste(type, " ", round(mean(vs<0)*100,1), " %"),
       xlab="Difference", col = colour)
  
  abline(v = 0, lty=2)
  
}

#stanplot(l2W_ranef, type = "hist",  pars="^b_")

# SPECIES number ------------------------------------------------------------------- 

sp_pois_f <- bf(SPEC_NO ~ 1 + TREAT + (1 + TREAT|GARDEN), family = poisson())
sp_prior <- get_prior(sp_pois_f, data = AllTestData)
#sp_prior$prior[1] <- "poisson()"

l3W_ranef <- brm(sp_pois_f, data=AllTestData,
           family = poisson(),
           control = list(adapt_delta = 0.999),
           prior = sp_prior,
           file = "spAllranef")


# plot(l3W_ranef)

pp = brms::pp_check(l3W_ranef)
pp + theme_bw()

sp_ranef_plot <- marginal_effects(l3W_ranef)

# Using predictions from the stanplots they produce the same results as predictions

## for comparison another plot using predictions of marginal_effect plot

# Plot template
ggplot(sp_ranef_plot$TREAT, aes(x = TREAT, y = estimate__)) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), position=position_dodge(), size=1, width=.5) +
  geom_jitter(data=AllTestData, aes(x=TREAT, y=SPEC_NO, group=GARDEN), 
              alpha=0.10, size=4,
              position = position_jitter(width = 0.07)) +
  geom_line(data=AllTestData, aes(x=TREAT, y=SPEC_NO, 
                                  group=GARDEN, 
                                  linetype = GARDEN), 
            size = 0.5, alpha=0.25) +
  geom_point(shape=21, size=4, fill='red') +
  xlab("") +
  ylab("Shannon's Index") +
  theme_bw () +
  theme(panel.grid = element_blank())

# Pairwise comparisons

sph <- stanplot(l3W_ranef, type = "hist", pars="^b_")

par(mfrow = c(3,2))

# Change the name of this

for(type in unique(sph$data$Parameter)){
  
  vs <- sph$data$value[sph$data$Parameter == type]
  
  colour = "lightblue"
  
  if (mean(vs<0) > 0.90) colour = "orange"
  if (mean(vs<0) > 0.95) colour = "red"
  
  hist(vs, breaks = 50,
       main = paste(type, " ", round(mean(vs<0)*100,1), " %"),
       xlab="Difference", col = colour)
  
  abline(v = 0, lty=2)
  
}

stanplot(l3W_ranef, type = "hist",  pars="^b_")


# Evenness -------------------------------------------------------------------

ev_norm_f <- bf(EVEN ~ 1 + TREAT + (1 + TREAT|GARDEN), family = Beta())
ev_prior <- get_prior(ev_norm_f, data = AllTestData)

l4W_ranef <- brm(ev_norm_f, data=AllTestData,
           family = Beta(),
           control = list(adapt_delta = 0.999),
           prior = ev_prior,
           file="eveAll")

pp = brms::pp_check(l4W_ranef)
pp + theme_bw()

ev_ranef_plot <- marginal_effects(l4W_ranef)

# Using predictions from the stanplots they produce the same results as predictions

## for comparison another plot using predictions of marginal_effect plot

# Plot template
ggplot(ev_ranef_plot$TREAT, aes(x = TREAT, y = estimate__)) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), position=position_dodge(), size=1, width=.5) +
  geom_jitter(data=AllTestData, aes(x=TREAT, y=EVEN, group=GARDEN), 
              alpha=0.10, size=4,
              position = position_jitter(width = 0.07)) +
  geom_line(data=AllTestData, aes(x=TREAT, y=EVEN, 
                                  group=GARDEN, 
                                  linetype = GARDEN), 
            size = 0.5, alpha=0.25) +
  geom_point(shape=21, size=4, fill='red') +
  xlab("") +
  ylab("Shannon's Index") +
  theme_bw () +
  theme(panel.grid = element_blank())

# Pairwise comparisons

sph <- stanplot(l4W_ranef, type = "hist", pars="^b_")

par(mfrow = c(3,2))

# Change the name of this

for(type in unique(sph$data$Parameter)){
  
  vs <- sph$data$value[sph$data$Parameter == type]
  
  colour = "lightblue"
  
  if (mean(vs<0) > 0.90) colour = "orange"
  if (mean(vs<0) > 0.95) colour = "red"
  
  hist(vs, breaks = 50,
       main = paste(type, " ", round(mean(vs<0)*100,1), " %"),
       xlab="Difference", col = colour)
  
  abline(v = 0, lty=2)
  
}


 
# Tree Dataset ---------------------------------------------------------------

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
