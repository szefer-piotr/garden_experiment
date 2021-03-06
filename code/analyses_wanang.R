# Examples
# https://www.sciencedirect.com/science/article/pii/S0095447017302310


# load data and packages --------------------------------------------------

# load data
AllTestData <- read.table("figs/stanplots/all_test_data.txt")

treat_levels = c("CONTROL", 
                 "FUNGICIDE",
                 "INSECTICIDE", 
                 "PREDATOR",
                 "WEEVIL25", 
                 "WEEVIL125")

AllTestData$TREAT = factor(AllTestData$TREAT,
                           levels = treat_levels)

TreeTestDataset <- read.table("figs/stanplots/tree_test_data.txt")
TreeTestDataset$TREATMENT = factor(TreeTestDataset$TREATMENT,
                           levels = treat_levels)

# Individual based biomass and trait data
inddat <- read.table("datasets/inddataWanang.txt")
#inddat$LIFE.FORM <-as.character(inddat$LIFE.FORM)
indtree <- inddat[inddat$LIFE.FORM %in% c("shrub","tree"), ]

indtree$TREAT = factor(indtree$TREAT, levels = treat_levels)
#indtree$SP_CODE = as.character(indtree$SP_CODE)


# Load packages
library(brms)
library(ggplot2)
library(bayesplot)
library(tidybayes)
library(dplyr)
library(ggridges)
library(nlme)
library(stringr)

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
                 control = list(adapt_delta = 0.8),
                 cores = 4, sample_prior = TRUE)

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

# [SUGESTION] Use more informative priors!

# Random intercept
# l1W <- brm(log(BIO)~TREAT + (1|GARDEN), 
#            data=AllTestData,
#            prior = priors_raninter2,
#            control = list(adapt_delta = 0.999))

# Specify priors for the random intercept and random slope model
priors_rsi2 <- c(set_prior("normal(4, 2)", class = "Intercept"),
                 set_prior("normal(0, 3)", class = "b"),
                 set_prior("normal(0, 4)", class = "sd"),
                 set_prior("normal(0, 3)", class = "sigma"),
                 set_prior("lkj(2)", class = "cor"))

## Whole community with random intercept and slope
l1W_ranef <- brm(log(BIO) ~ TREAT + (1+TREAT|GARDEN), 
                 data=AllTestData,
                 prior = priors_rsi2,
                 control = list(adapt_delta = 0.999),
                 file="liW_ranef")

# Comparisons with model without random slope.
# compare_ic(waic(l1W), waic(l1W_ranef), ic = "waic")

# Model diagnostics
pp = brms::pp_check(l1W_ranef)
pp + theme_bw()


# plots

# Old plot
# AllTestData %>% 
#   tidybayes::add_predicted_draws(l1W_ranef, n = 1000) %>%
#   ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=TREAT)) + 
#   geom_point(aes(x = log(BIO), y = TREAT), size = 3, alpha = 0.1) +
#   geom_density_ridges(alpha = 0.5)

# Pairwise comparisons

# Included in the suplementary and determining "significance"
sph <- stanplot(l1W_ranef, type = "hist", pars="^b_")

png("figs/bio_all_diff.png")

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

dev.off()

# Colours based on posterior distributions of differences between
# control and a treatment
colour = as.character(AllTestData$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "grey50"
colour[colour == "PREDATOR"] <- "orange"
colour[colour == "WEEVIL25"] <- "grey50"
colour[colour == "WEEVIL125"] <- "red"

AllTestData$colour = colour

png("figs/bio_bayes_hist.png")
AllTestData %>% 
  tidybayes::add_predicted_draws(l1W_ranef, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = log(BIO), y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()
dev.off()

#### Hold the plot
BD <- AllTestData %>% 
  tidybayes::add_predicted_draws(l1W_ranef, n = 1000)

BDp <- ggplot(BD,  aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = log(BIO), y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()
####


# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html

garden_panel_plot <- function(x) x %>%
  ggplot(aes(x = TREAT, y = log(BIO))) +
  geom_point(aes(y = .prediction), alpha = 0.1, position = position_jitter(width = 0.1)) +
  geom_point(colour = "red", size = 3) +
  facet_wrap(~GARDEN) +
  coord_flip() +
  theme_minimal()

# predictions for the AVERAGE block (supplementary?)
AllTestData %>%
  tidybayes::add_predicted_draws(l1W_ranef, n = 200, re_formula = NULL) %>%
  garden_panel_plot



# 1b. Biomass - random species --------------------------------------------

# https://rpsychologist.com/r-guide-longitudinal-lme-lmer


# Specify the equation
# UNDER CONSTRUCTION - DON'T EVEN KNOW IF THIS IS OK

bio_ln_ind_rbsp <- bf(WEIGHT ~ 1 + TREAT + (1 + TREAT|BLOCK) + (1 + TREAT|SP_CODE),
                   family = "lognormal")

bio_ln_ind_rintsp <- bf(WEIGHT ~ 1 + TREAT + (1 + TREAT|BLOCK) + (1|SP_CODE))

bio_ln_ind_rb <- bf(WEIGHT ~ 1 + TREAT + (1 + TREAT|BLOCK),
                    family = "lognormal")

# nlme estimates
blirbsp <- lme

## Set some priors 
priors <- c(set_prior("normal(4, 2)", class = "Intercept"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 4)", class = "sd"),
            set_prior("normal(0, 3)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

lnbio_random_block <- brm(bio_ln_ind_rb,
                    data=indtree,
                    control = list(adapt_delta = .8))

lnbio_random_block_int_sp <- brm(bio_ln_ind_rintsp,
                                  data=indtree,
                                  control = list(adapt_delta = .8))

lnbio_random_block_species <- brm(bio_ln_ind_rbsp,
                    data=indtree,
                    control = list(adapt_delta = .8))

compare_ic(waic(lnbio_random_block), 
           waic(lnbio_random_block_species), 
           waic(lnbio_random_block_int_sp),
           ic = "waic")

sph <- stanplot(lnbio_random_block, type = "hist", pars="^b_")

png("figs/bio_ind_tree_diff.png")
X11(width=10, height=10)
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
dev.off()

colour = as.character(indtree$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "grey50"
colour[colour == "PREDATOR"] <- "grey50"
colour[colour == "WEEVIL25"] <- "red"
colour[colour == "WEEVIL125"] <- "red"
indtree$colour = colour

#png("figs/bio_ind_bayes_hist.png")
indtree %>% 
  tidybayes::add_predicted_draws(lnbio_random_block_int_sp, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "red","orange")) +
  #geom_point(aes(x = log(WEIGHT), y = TREAT), size = 3, alpha = 0.1,
  #           color = "grey50") + 
  theme_bw()
#dev.off()

# 1c. leaf damage ----------------------

# Leaf damage logit values
ld_gaus_ind <- bf(HERB ~ 1 + TREAT + (1 + TREAT|BLOCK) + (1 + TREAT|SP_CODE))

ld_prior <- c(set_prior("normal(-10, 4)", class = "Intercept"),
              set_prior("normal(0, 3)", class = "b"),
              set_prior("normal(0, 4)", class = "sd"),
              set_prior("normal(0, 3)", class = "sigma"),
              set_prior("lkj(2)", class = "cor"))

ld_rand_ind <- brm(ld_gaus_ind,
                   data=indtree,
                   control = list(adapt_delta = .8),
                   prior = ld_prior)


# Transform them back to proportions
library(boot)
indtree$HERB <- inv.logit(indtree$HERB)
indtree$HERB[indtree$HERB == 0] <- 1e-07

ld_beta_ind <- bf(HERB ~ 1 + TREAT + (1 + TREAT|BLOCK) + (1 + TREAT|SP_CODE),
                   family = Beta())

ld_rand_ind <- brm(ld_beta_ind,
                    data=indtree,
                    control = list(adapt_delta = .8))

stanplot(ld_rand_ind)

pp = brms::pp_check(ld_rand_ind)
pp + theme_bw()

sph <- stanplot(ld_rand_ind, type = "hist", pars="^b_")

#png("figs/ld_tree_diff_ind.png")

X11(10,10)
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
#dev.off()

indtree %>% 
  tidybayes::add_predicted_draws(ld_rand_ind, n = 1000) %>%
  ggplot(aes(x = logit(.prediction), y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "red","orange")) +
  #geom_point(aes(x = log(WEIGHT), y = TREAT), size = 3, alpha = 0.1,
  #           color = "grey50") + 
  theme_bw()

# 2. SHANNON ------------------------------------------------------------------

sw_norm_f <- bf(SW ~ 1 + TREAT + (1 + TREAT|GARDEN), family = gaussian())

mean(AllTestData$SW)

sw_prior <- c(set_prior("normal(1, 3)", class = "Intercept"),
            set_prior("normal(0, 3)", class = "b"),
            set_prior("normal(0, 4)", class = "sd"),
            set_prior("normal(0, 3)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

# How to make it better? They shouldn't be smaller tan zero!
l2W_ranef <- brm(SW ~ TREAT + (1+TREAT|GARDEN),
           data=AllTestData,
           control = list(adapt_delta = 0.999),
           prior = sw_prior)


#plot(l2W_ranef)

pp = brms::pp_check(l2W_ranef)
pp + theme_bw()


# Pairwise comparisons

sph <- stanplot(l2W_ranef, type = "hist", pars="^b_")

png("figs/sw_all_diff.png")
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
dev.off()

colour = as.character(AllTestData$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "red"
colour[colour == "PREDATOR"] <- "grey50"
colour[colour == "WEEVIL25"] <- "orange"
colour[colour == "WEEVIL125"] <- "grey50"

AllTestData$colour = colour

png("figs/sw_bayes_hist.png")
AllTestData %>% 
  tidybayes::add_predicted_draws(l2W_ranef, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = SW, y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()
dev.off() 

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

#stanplot(l2W_ranef, type = "hist",  pars="^b_")

# Checkong model predictions... shouldn't predict values smaller than zero!
garden_panel_plot <- function(x) x %>%
  ggplot(aes(x = TREAT, y = SW)) +
  geom_point(aes(y = .prediction), alpha = 0.1, position = position_jitter(width = 0.1)) +
  geom_point(colour = "red", size = 3) +
  facet_wrap(~GARDEN) +
  coord_flip() +
  theme_minimal()

# predictions for the AVERAGE block (supplementary?)
AllTestData %>%
  tidybayes::add_predicted_draws(l2W_ranef, n = 200, re_formula = NULL) %>%
  garden_panel_plot
# end of stuff


# 3. SPECIES number ------------------------------------------------------------------- 

sp_negb_f <- bf(SPEC_NO ~ 1 + TREAT + (1 + TREAT|GARDEN), family = "negbinomial")
sp_pois_f <- bf(SPEC_NO ~ 1 + TREAT + (1 + TREAT|GARDEN), family = "poisson")

# Set some priors

# I DONT KNOW HOW TO SPECIFY THEM!!!
# sp_negbin_priors <- c(
#   set_prior("poisson(15)", class = "b"),
#   set_prior("lkj(1)", class = "cor"),
#   set_prior("poisson(15)", class = "Intercept"),
#   set_prior("normal(0,1)", class = "sd")
# )

# Prior predictive check!
# sample from the prior distribution -> # sample_prior = "only"

l3W_ranef_negbin <- brm(sp_negb_f, data=AllTestData,
                        family = "negbinomial",
                        control = list(adapt_delta = 0.8))

l3W_ranef_pois <- brm(sp_pois_f, data=AllTestData,
                      family = "poisson",
                      control = list(adapt_delta = 0.8))

# I guess that Poisson stays?
compare_ic(waic(l3W_ranef_pois), waic(l3W_ranef_negbin), ic = "waic")

colour = as.character(AllTestData$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "orange"
colour[colour == "INSECTICIDE"] <- "grey50"
colour[colour == "PREDATOR"] <- "grey50"
colour[colour == "WEEVIL25"] <- "grey50"
colour[colour == "WEEVIL125"] <- "grey50"

AllTestData$colour = colour

png("figs/rich_bayes_hist.png")
AllTestData %>% 
  tidybayes::add_predicted_draws(l3W_ranef_pois, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = SPEC_NO, y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()
dev.off()

png("figs/rich_all_diff.png")
par(mfrow = c(3,2))

sph <- stanplot(l3W_ranef_pois, type = "hist", pars="^b_")
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
dev.off()

# Using predictions from the stanplots they produce the same results as predictions
## for comparison another plot using predictions of marginal_effect plot


# 4. Evenness -------------------------------------------------------------------

ev_norm_f <- bf(EVEN ~ 1 + TREAT + (1 + TREAT|GARDEN), family = Beta())

# ??? 
#ev_prior <- c()

l4W_ranef <- brm(ev_norm_f, data=AllTestData,
           family = "Beta",
           control = list(adapt_delta = 0.8), file = "l4W_ranef")


pp = brms::pp_check(l4W_ranef)
pp + theme_bw()


sph <- stanplot(l4W_ranef, type = "hist", pars="^b_")
png("figs/eve_all_diff.png")
#X11(10,10)
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
dev.off()


colour = as.character(AllTestData$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "orange"
colour[colour == "PREDATOR"] <- "grey50"
colour[colour == "WEEVIL25"] <- "orange"
colour[colour == "WEEVIL125"] <- "grey50"

AllTestData$colour = colour

png("figs/eve_bayes_hist.png")
AllTestData %>% 
  tidybayes::add_predicted_draws(l4W_ranef, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = EVEN, y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()
dev.off()



# 5. Tree Biomass ---------------------------------------------------------------

TreeTestDataset$TREATMENT


# Specify priors
bio_norm_tree <- bf(logBio ~ 1 + TREATMENT + (1 + TREATMENT|GARDEN), family = gaussian())

# Model
l1Wt <- brm(bio_norm_tree, 
            data=TreeTestDataset,
            control = list(adapt_delta = 0.99),
            file="bioTre")

pp = brms::pp_check(l1Wt)
pp + theme_bw()


sph <- stanplot(l1Wt, type = "hist", pars="^b_")
png("figs/bio_tree_diff.png")
#X11(10,10)
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
dev.off()


colour = as.character(TreeTestDataset$TREAT)
colour[colour == "CONTROL"] <- "grey50"
colour[colour == "FUNGICIDE"] <- "grey50"
colour[colour == "INSECTICIDE"] <- "grey50"
colour[colour == "PREDATOR"] <- "grey50"
colour[colour == "WEEVIL25"] <- "grey50"
colour[colour == "WEEVIL125"] <- "red"

TreeTestDataset$colour = colour

png("figs/bio_bayes_tree_hist.png")
TreeTestDataset %>% 
  tidybayes::add_predicted_draws(l1Wt, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREATMENT, group = TREATMENT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "red","orange")) +
  geom_point(aes(x = logBio, y = TREATMENT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()
dev.off()

# 6. Tree Diversity --------------------------------------------------------------

# l2Wtbio <- brm(DivTree~TREATMENT + (1 + TREATMENT|GARDEN) + logBio, data=TreeTestDataset,
#             control = list(adapt_delta = 0.99))

lognormDiv = bf(DivTree~TREATMENT + (1 + TREATMENT|GARDEN), family = lognormal())

l2Wt_ran <- brm(lognormDiv, data=TreeTestDataset,
               control = list(adapt_delta = 0.99))

# lognorm_int_Div = bf(DivTree~TREATMENT + (1|GARDEN), family = "lognormal")
# 
# l2Wt_int <- brm(lognorm_int_Div, data=TreeTestDataset,
#                control = list(adapt_delta = 0.99))


# l2Wt_invSimp <- brm(invSimp~TREATMENT + (1 + TREATMENT|GARDEN), data=TreeTestDataset,
#             control = list(adapt_delta = 0.8))

compare_ic(waic(l2Wt_ran), waic(l2Wt_int), ic = "waic")

pp <- brms::pp_check(l2Wt_ran)
pp + theme_bw()

TreeTestDataset %>% 
  tidybayes::add_predicted_draws(l2Wt_ran, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREATMENT, group = TREATMENT)) + 
  geom_density_ridges(alpha = 0.5) + xlim(-1,3) +
  scale_fill_manual(values = c("grey", "red","orange")) +
  geom_point(aes(x = DivTree, y = TREATMENT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()


# Extreme predictions
ep = TreeTestDataset %>% 
  tidybayes::add_predicted_draws(l2Wt_ran, n = 1000)

# Hey, I think you better use Hill numbers (so exponential Shannon and inverse Simpson) for diversity calculations and comparisons. Cf. Jost, L. (2006). Entropy and diversity. Oikos, 113(2), 363–375.
# or: Chao, A., Chiu, C.-H., & Jost, L. (2014). Unifying Species Diversity, Phylogenetic Diversity, Functional Diversity, and Related Similarity and Differentiation Measures Through Hill Numbers. 
# Annual Review of Ecology, Evolution, and Systematics, 45, 297–324. doi:10.1146/annurev-ecolsys-120213-091540

# There are good and bad diversity measurement concepts and the distinction between them is determined by properties, contexts and individual judgements.
# I suggest you try to asses the "true" species diversity, using species richness and the effective number of species, sensu Jost (2006).
# Jost, L. 2006. Entropy and diversity. - Oikos 113: 363–374

hist(log(ep$.prediction))

# 7. Richness ---------------------------------
rich_bf = bf(species~TREATMENT+(1+TREATMENT|GARDEN), family = "poisson")

l3Wt <- brm(rich_bf, 
            data=TreeTestDataset,
            control = list(adapt_delta = 0.8))

TreeTestDataset %>%
  tidybayes::add_predicted_draws(l3Wt, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREATMENT, group = TREATMENT)) + 
  geom_density_ridges(alpha = 0.5) +
  scale_fill_manual(values = c("grey", "red","orange")) +
  geom_point(aes(x = species, y = TREATMENT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()

sph <- stanplot(l3Wt, type = "hist", pars="^b_")
#png("figs/bio_tree_diff.png")
X11(10,10)
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

dev.off()
  
# 8. Evenness ---------------
# l4Wt <- brm(Evenness~TREATMENT + (1|GARDEN), data=TreeTestDataset,
#              family = Beta(),
#              control = list(adapt_delta = 0.99),
#              file = "eveTre")

# 9. Stem number -----------------

l5Wt <- brm(stems~TREATMENT + (1+TREATMENT|GARDEN), 
            data=TreeTestDataset,
            family=poisson(),
            control = list(adapt_delta = 0.8))


l5Wt_nb <- brm(stems~TREATMENT + (1+TREATMENT|GARDEN), 
            data=TreeTestDataset,
            family= "negbinomial",
            control = list(adapt_delta = 0.8))

compare_ic(waic(l5Wt), waic(l5Wt_nb), ic = "waic")

hist(TreeTestDataset$stems)


TreeTestDataset %>%
  tidybayes::add_predicted_draws(l5Wt, n = 1000) %>%
  ggplot(aes(x = .prediction, y = TREATMENT, group = TREATMENT)) + 
  geom_density_ridges(alpha = 0.5) + xlim(c(0,100)) +
  scale_fill_manual(values = c("grey", "red","orange")) +
  geom_point(aes(x = stems, y = TREATMENT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()


sph <- stanplot(l5Wt, type = "hist", pars="^b_")
#png("figs/bio_tree_diff.png")
X11(10,10)
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


# 9.b stems per individual --------

l5Wt_ind <- brm(NO_STEMS~TREAT + (1+TREAT|BLOCK) + (1+TREAT|SP_CODE), 
            data=indtree,
            family=poisson(),
            control = list(adapt_delta = 0.8))



indtree %>%
  tidybayes::add_predicted_draws(l5Wt_ind, n = 1000, allow_new_levels = TRUE) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT)) + 
  geom_density_ridges(alpha = 0.5) + xlim(c(0,100)) +
  scale_fill_manual(values = c("grey", "red","orange")) +
  geom_point(aes(x = NO_STEMS, y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()

# 10. SLA --------------

sla <- brm(log(SLA) ~ TREAT + (1+TREAT|BLOCK) + (1+TREAT|SP_CODE),
           data = inddat,
           control = list(adapt_delta = 0.8))

sph <- stanplot(sla, type = "hist", pars="^b_")
X11(10,10)
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



inddat %>% 
  tidybayes::add_predicted_draws(sla, n = 1000, allow_new_levels = TRUE) %>%
  ggplot(aes(x = .prediction, y = TREAT, group = TREAT)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "red","orange")) +
  #geom_point(aes(x = log(WEIGHT), y = TREAT), size = 3, alpha = 0.1,
  #           color = "grey50") + 
  theme_bw()

# 11. LDMC -----------------

# Plots in a panel ------------------

#### Hold the plot
BD <- AllTestData %>% 
  tidybayes::add_predicted_draws(l1W_ranef, n = 1000)

BDp <- ggplot(BD,  aes(x = .prediction, y = TREAT, group = TREAT, fill=colour)) + 
  geom_density_ridges(alpha = 0.5) + 
  scale_fill_manual(values = c("grey", "orange", "red")) +
  geom_point(aes(x = log(BIO), y = TREAT), size = 3, alpha = 0.1,
             color = "grey50") + 
  theme_bw()
####