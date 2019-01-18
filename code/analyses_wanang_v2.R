# 0. Load data and packages ----

library(dplyr)
library(lme4)
library(lmerTest)
library(Hmisc) #for the stat_summary plots

# Main dataset
main <- read.table("datasets/wng_main_clean.txt", header = T)

# Tree dataset
tree <- main %>%
  filter(LIFE.FORM %in% c("tree", "shrub"))

# Descriptive statistics
test <- read.table("datasets/wng_main_test.txt")
tree_test <- read.table("datasets/wng_tree_test.txt")

# 0.b Remove insecticide ----
test <- test[rownames(test) != "WG1P6", ]
tree_test <- tree_test[rownames(tree_test) != "WG1P6", ]
ggplot(test, aes(y=log(BIO), x=TREAT)) + geom_point()


# 1. BIOMASS models -----

# Cumulative biomass
bio_rbl <- lmer(log(BIO) ~ TREAT + (1|GARDEN),
                data = test) 

summary(bio_rbl)

# Richness
spn_rbl <- lmer(SPEC_NO ~ TREAT + (1|GARDEN),
                data = test)

summary(spn_rbl)

# Diversity
div_rbl <- lmer(SW ~ TREAT + (1|GARDEN),
                data = test)

summary(div_rbl)

# Number of stems
abu_rbl <- glmer(stems ~ TREATMENT + (1|GARDEN),
                data = tree_test, family = "poisson")

abu_nb_rbl <- glmer.nb(stems ~ TREATMENT + (1|GARDEN),
                 data = tree_test)


summary(abu_rbl)
summary(abu_nb_rbl)

library(blmeco)
# Here is the test for overdispersion and with Poisson there is an overdispersion
dispersion_glmer(abu_rbl) # should not exceed 1.4
dispersion_glmer(abu_nb_rbl)

AIC(abu_rbl, abu_nb_rbl)

# Optional likelihood ratio test
# #compute a model where the effect of status is estimated
# unrestricted_fit <- lmer(
#   formula = SPEC_NO ~ (1|GARDEN) + TREAT,
#   REML = F,
#   data = test
# )
# 
# restricted_fit = lmer(
#   formula = SPEC_NO ~ (1|GARDEN),
#   REML = F,
#   data = test#because we want to compare models on likelihood
# )
# 
# #compute the AIC-corrected log-base-2 likelihood ratio (a.k.a. "bits" of evidence)
# (AIC(restricted_fit)-AIC(unrestricted_fit))*log2(exp(1))

# Generate plots! -----
# Biomass for all, shannon's diverstiy, number of stems
panel_plot_staked_data <- rbind(setNames(test[,c("TREAT","BIO")],c("treat","val")),
                                setNames(test[,c("TREAT","SW")],c("treat","val")),
                                setNames(tree_test[,c("TREATMENT","stems")], c("treat","val")))

# Refactor
panel_plot_staked_data$treat <- factor(panel_plot_staked_data$treat,
                                       labels = c("C", "F", "I", "P", "H2", "H1"))
panel_plot_staked_data$treat  <- ordered(panel_plot_staked_data$treat, 
                 levels = c("C", "F", "I", "P", "H1", "H2"))

panel_plot_staked_data$type <- rep(c("Biomass", "Diversity", "Abundance"), each = dim(test)[1])

panel_plot_staked_data$type <- ordered(panel_plot_staked_data$type, 
                                         levels = c("Biomass", "Diversity", "Abundance"))

# >>> FIG1 -----
fig1 <- ggplot(panel_plot_staked_data, aes(x=treat, y=val)) + 
  facet_wrap(~type, scales = "free") + 
  geom_point(col = "grey80", alpha = 0.5) + 
  theme_bw() +
  theme(axis.text.x=element_text(angle=0, hjust=0.5),
        strip.text = element_text(size=10),
        legend.justification=c(0.5,0.5), 
        legend.position="bottom")+xlab("")+ylab("")

colors = color=c("black","grey50","grey50","grey50","grey50","red",
                 "black","grey50","red","grey50","grey50","grey50",
                 "black","red","grey50","grey50","grey50","red")

fig1 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color= colors, width=0.2, lwd=0.7) +
  stat_summary(fun.y=mean, geom="point", color=colors, cex = 3)
  #scale_y_continuous(breaks = round(seq(0, 300, by = 50),1))
