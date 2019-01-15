main <- read.table("datasets/wng_ind_dbh.txt", header = T)

library(dplyr)

main %>%
  filter(LIFE.FORM %in% c("tree", "shrub") & NO_STEMS != "NA") %>%
  group_by(SP_CODE)
