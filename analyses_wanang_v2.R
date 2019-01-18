# 0. Load data and packages ----

library(dplyr)
library(lme4)
library(lmerTest)

# Main dataset
main <- read.table("datasets/wng_main.txt", header = T)

# Tree dataset
tree <- main %>%
  filter(LIFE.FORM %in% c("tree", "shrub"))

# Descriptive statistics
test <- read.table("figs/stanplots/all_test_data.txt")
  
# 1. BIOMASS models -----

# Singular fit!

# Per species!
# sp_bio_rbl <- lmer(log(WEIGHT) ~ TREAT + (1|BLOCK),
#                    data = main)
#summary(sp_bio_rbl)

# Cumulative biomass
bio_rbl <- lmer(log(BIO) ~ TREAT + (1|GARDEN),
                data = test) 

summary(bio_rbl)

# Diversity
div_rbl <- lmer(SW ~ TREAT + (1|GARDEN),
                data = test)

summary(div_rbl)
