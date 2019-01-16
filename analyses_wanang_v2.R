# 0. Load data and packages ----

library(dplyr)
library()

# Main dataset
main <- read.table("datasets/wng_main.txt", header = T)
# Tree dataset
tree <- main %>%
  filter(LIFE.FORM %in% c("tree", "shrub"))

# Descriptive statistics
test <- read.table("figs/stanplots/all_test_data.txt")
  
# 1. BIOMASS models -----

