

# Examples
# https://www.sciencedirect.com/science/article/pii/S0095447017302310


# load data and packages --------------------------------------------------

# load data
AllTestData <- read.table("figs/stanplots/all_test_data.txt")

AllTestData$TREAT = factor(AllTestData$TREAT,
                           levels = c("CONTROL", 
                                      "FUNGICIDE",
                                      "INSECTICIDE", 
                                      "PREDATOR",
                                      "WEEVIL25", 
                                      "WEEVIL125"))


TreeTestDataset <- read.table("figs/stanplots/tree_test_data.txt")

# Load package 
library(brms)
library(ggplot2)
library(bayesplot)
library(tidybayes)
library(dplyr)
library(ggforce)

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
marginal_effects(l1W_ranef_lnorm)

# recoding factors --------------------------------------------------------

#the model above still has a hard time with divergent iterations. I am curious
#to try how it works if the variables are coded by hand to develop contrasts
library(stringr)

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


## 1. Biomass -----------------------------------------------------

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

# model diagnostics

pp = brms::pp_check(l1W_ranef)
pp + theme_bw()


# plots
library(ggridges)

AllTestData %>% 
  tidybayes::add_predicted_draws(l1W_ranef, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=TREAT)) + 
  geom_point(aes(x = log(BIO), y = TREAT), size = 3, alpha = 0.1) +
  geom_density_ridges(alpha = 0.5)

colour = as.character(AllTestData$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "grey50"
colour[colour == "PREDATOR"] <- "orange"
colour[colour == "WEEVIL25"] <- "grey50"
colour[colour == "WEEVIL125"] <- "red"

AllTestData$colour = colour

AllTestData %>% 
  tidybayes::add_predicted_draws(l1W_ranef, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = log(BIO), y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()


# Pairwise comparisons

# Included in the suplementary and determining "significance"
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

# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

garden_panel_plot <- function(x) x %>%
  ggplot(aes(x = TREAT, y = log(BIO))) +
  geom_point(aes(y = .prediction), alpha = 0.1, position = position_jitter(width = 0.1)) +
  geom_point(colour = "red", size = 3) +
  facet_wrap(~GARDEN) +
  coord_flip() +
  theme_minimal()


# predictions for the AVERAGE block
AllTestData %>%
  tidybayes::add_predicted_draws(l1W_ranef, n = 200, re_formula = NULL) %>%
  garden_panel_plot

# predictions for the AVERAGE block
new_plot_data <- AllTestData %>%
  modelr::data_grid(TREAT, GARDEN = "PS1") %>%
  tidybayes::add_predicted_draws(l1W_ranef, n = 200, re_formula = NA)


AllTestData %>%
  ggplot(aes(x = TREAT, y = log(BIO))) +
  geom_point(aes(y = .prediction), alpha = 0.1, position = position_jitter(width = 0.1),
             new_plot_data) +
  geom_point(colour = "green", size = 3) +
  facet_wrap(~GARDEN) +
  coord_flip() +
  theme_minimal()


# Histogram for the differences

# 1b. Biomass - random species --------------------------------------------

inddat <- read.table("datasets/inddataWanang.txt")
inddat$LIFE.FORM <-as.character(inddat$LIFE.FORM)

indtree <- inddat[inddat$LIFE.FORM %in% c("shrub","tree"), ]

# https://rpsychologist.com/r-guide-longitudinal-lme-lmer
# Specify the equation

bio_lnorm_ind <- bf(WEIGHT ~ 1 + TREAT + (1 + SP_CODE + TREAT|BLOCK), family = lognormal())

bio_rand_ind <- brm(bio_lnorm_ind,
                    data=indtree,
                    control = list(adapt_delta = .8))



# 2. SHANNON ------------------------------------------------------------------

sw_norm_f <- bf(SW ~ 1 + TREAT + (1 + TREAT|GARDEN), family = gaussian())

mean(AllTestData$SW)

sw_prior <- c(set_prior("normal(0, 3)", class = "Intercept"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 4)", class = "sd"),
            set_prior("normal(0, 3)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

# How to get better 
l2W_ranef <- brm(SW ~ TREAT + (1+TREAT|GARDEN),
           data=AllTestData,
           control = list(adapt_delta = 0.999),
           prior = sw_prior,
           file = "swAllranef")


#plot(l2W_ranef)

pp = brms::pp_check(l2W_ranef)
pp + theme_bw()

## for comparison another plot using predictions of marginal_effect plot
# Plot template

colour = as.character(AllTestData$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "red"
colour[colour == "PREDATOR"] <- "grey50"
colour[colour == "WEEVIL25"] <- "orange"
colour[colour == "WEEVIL125"] <- "grey50"

AllTestData$colour = colour

AllTestData %>% 
  tidybayes::add_predicted_draws(l2W_ranef, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = SW, y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()
  

# sw_ranef_plot <- marginal_effects(l2W_ranef)    # second model seems to better represent the data

# ggplot(sw_ranef_plot$TREAT, aes(x = TREAT, y = estimate__)) +
#   geom_errorbar(aes(ymin=lower__, ymax=upper__), position=position_dodge(), size=1, width=.5) +
#   geom_jitter(data=AllTestData, aes(x=TREAT, y=SW, group=GARDEN), 
#               alpha=0.10, size=4,
#               position = position_jitter(width = 0.07)) +
#   geom_line(data=AllTestData, aes(x=TREAT, y=SW, 
#                                   group=GARDEN, 
#                                   linetype = GARDEN), 
#             size = 0.5, alpha=0.25) +
#   geom_point(shape=21, size=4, fill='red') +
#   xlab("") +
#   ylab("Shannon's Index") +
#   theme_bw () +
#   theme(panel.grid = element_blank())

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

# 3. SPECIES number ------------------------------------------------------------------- 

sp_pois_f <- bf(SPEC_NO ~ 1 + TREAT + (1 + TREAT|GARDEN), family = "negbinomial")

# Set some priors

# <<<<< I DONT KNOW HOW TO SPECIFY THEM!!! >>>>

# sp_no_priors <- c(
#   set_prior("poisson(15)", class = "b"),
#   set_prior("lkj(1)", class = "cor"),
#   set_prior("poisson(15)", class = "Intercept"),
#   set_prior("normal(0,1)", class = "sd")
# )

# Prior predictive check!
# sample from the prior distribution
l3W_ranef_prior <- brm(sp_pois_f, data=AllTestData,
           family = "negbinomial",
           control = list(adapt_delta = 0.8),
           prior = sp_no_priors, sample_prior = "only")

l3W_ranef_no_prior <- brm(sp_pois_f, data=AllTestData,
                       family = "negbinomial",
                       control = list(adapt_delta = 0.8),
                       file = "l3W_ranef_no_prior")


l3W_ranef_no_prior %>% summary


####  3b. Plot ----------------------------------------------------

pp = brms::pp_check(l3W_ranef_no_prior)
pp + theme_bw()

## for comparison another plot using predictions of marginal_effect plot


# Pairwise comparisons

sph <- stanplot(l3W_ranef_no_prior, type = "hist", pars="^b_")

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



# Plot template

colour = as.character(AllTestData$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "grey50"
colour[colour == "PREDATOR"] <- "grey50"
colour[colour == "WEEVIL25"] <- "grey50"
colour[colour == "WEEVIL125"] <- "grey50"

AllTestData$colour = colour

AllTestData %>% 
  tidybayes::add_predicted_draws(l3W_ranef_no_prior, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = SPEC_NO, y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()



sph <- stanplot(l3W_ranef_prior, type = "hist", pars="^b_")
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


# visualize the prior predictive distribution:
AllTestData %>% 
  tidybayes::add_predicted_draws(l3W_ranef_no_prior) %>% #glimpse %>% 
  filter(.draw ==3) %>% 
  ggplot(aes(x = TREAT, y = .prediction)) + geom_point()

# sample from the prior distribution
l3W_ranef <- brm(sp_pois_f, data=AllTestData,
                       family = poisson(),
                       control = list(adapt_delta = 0.8),
                       prior = sp_no_priors, sample_prior = "yes")


# 1. Visualize fluffy plots
AllTestData %>% 
  tidybayes::add_predicted_draws(l3W_ranef_no_prior, n = 500) %>% glimpse %>% 
  ggplot(aes(x = TREAT, y = SPEC_NO)) + 
  geom_point(aes(y = .prediction), alpha = 0.1, position = position_jitter(width = 0.1)) + 
  geom_point(colour = "green", size = 3) + 
  facet_wrap(~GARDEN) + 
  coord_flip() + 
  theme_minimal()

pp_check(l3W_ranef, fun = "stat_grouped")
pp_check

stanplot(l3W_ranef, type = "dens")

launch_shinystan(l3W_ranef)

# 2. Change the family


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


# 4. Evenness -------------------------------------------------------------------

ev_norm_f <- bf(EVEN ~ 1 + TREAT + (1 + TREAT|GARDEN), family = Beta())

# ??? 
#ev_prior <- c()

l4W_ranef <- brm(ev_norm_f, data=AllTestData,
           family = "Beta",
           control = list(adapt_delta = 0.8), file = "l4W_ranef")


# 4b. Plots --------------------

pp = brms::pp_check(l4W_ranef)
pp + theme_bw()


colour = as.character(AllTestData$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "orange"
colour[colour == "PREDATOR"] <- "grey50"
colour[colour == "WEEVIL25"] <- "orange"
colour[colour == "WEEVIL125"] <- "grey50"

AllTestData$colour = colour

AllTestData %>% 
  tidybayes::add_predicted_draws(l4W_ranef, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = EVEN, y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()



# Using predictions from the stanplots they produce the same results as predictions

## for comparison another plot using predictions of marginal_effect plot

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


 

# Tree Biomass ---------------------------------------------------------------

# Specify priors
bio_norm_tree <- bf(logBio ~ 1 + TREATMENT + (1 + TREATMENT|GARDEN), family = gaussian())
bio_tree_prior <- get_prior(bio_norm_tree, data = TreeTestDataset)

# Model
l1Wt <- brm(bio_norm_tree, 
            data=TreeTestDataset,
            control = list(adapt_delta = 0.99),
            prior = bio_tree_prior,
            file="bioTre")


bio_ranef_tree <- marginal_effects(l1Wt)

## for comparison another plot using predictions of marginal_effect plot
# Plot template

p <- ggplot(bio_ranef_tree$TREAT, aes(x = TREATMENT, y = estimate__)) +
  geom_errorbar(aes(ymin=lower__, ymax=upper__), position=position_dodge(), size=1, width=.5) +
  geom_jitter(data=TreeTestDataset, aes(x=TREATMENT, y=logBio, group=GARDEN), 
              alpha=0.10, size=4,
              position = position_jitter(width = 0.07)) +
  geom_line(data=TreeTestDataset, aes(x=TREATMENT, y=logBio, 
                                  group=GARDEN, 
                                  linetype = GARDEN), 
            size = 0.5, alpha=0.25) +
  geom_point(shape=21, size=4, fill='red') +
  xlab("") +
  ylab("LN(BIOMASS)") +
  theme_bw () +
  theme(panel.grid = element_blank())

p

# plotMyResults(l1Wt)

# Tree Diversity --------------------------------------------------------------

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





#### Under construction!
plotMyResults <- function(l2W_ranef){
  
  pp = brms::pp_check(l2W_ranef)
  pp + theme_bw()
  
  sw_ranef_plot <- marginal_effects(l2W_ranef)    # second model seems to better represent the data
  
  ## for comparison another plot using predictions of marginal_effect plot
  # Plot template
  
  p <- ggplot(sw_ranef_plot$TREAT, aes(x = TREAT, y = estimate__)) +
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
  
  p
  
}
