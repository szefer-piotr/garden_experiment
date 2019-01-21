# 0. Load data and packages ----

library(dplyr)
library(lme4)
library(lmerTest)
library(Hmisc) #for the stat_summary plots
library(brms)


contingencyTable2 <- function(dataset, ROW, COL, VALUE){
  # Get rid of the empty factors
  dataset[, colnames(dataset) == ROW] <- as.character(dataset[, colnames(dataset) == ROW])
  dataset[, colnames(dataset) == COL] <- as.character(dataset[, colnames(dataset) == COL])
  # Make a table, get rid of the empty rows and columns
  plants <- table(dataset[, colnames(dataset) == ROW], dataset[, colnames(dataset) == COL])
  plants <- plants[rowSums(plants) != 0, colSums(plants) != 0]
  # See where to insert values
  allSpecCodes <- colnames(plants)
  allPlotCodes <- rownames(plants)
  entries <- which(plants != 0, arr.ind = TRUE)
  # Loop through the entries and insert values
  for (entry in 1:dim(entries)[1]){
    plot <- entries[entry,1]
    plant <- entries[entry,2]
    specCode <- allSpecCodes[plant]
    plotCode <- allPlotCodes[plot]
    #res <- dataset[dataset$ROW == plotCode & dataset$COL == specCode,VALUE]
    res <- dataset[dataset[,ROW] == plotCode & dataset[,COL] == specCode,VALUE]
    # print(sum(res))
    plants[plot,plant] <- sum(res, na.rm = TRUE)
  }
  plants[is.na(plants)] <- 0
  
  # Change the table to a matrix or data.frame
  
  mat_a <- matrix(0, nrow = dim(plants)[1], ncol = dim(plants)[2])
  colnames(mat_a) <- colnames(plants)
  rownames(mat_a) <- rownames(plants)
  for (row in 1:dim(plants)[1]){
    for (col in 1:dim(plants)[2])
      mat_a[row,col] <- plants[row,col]
  }
  
  #return(plants)
  return(mat_a)
}

# Main dataset
main <- read.table("datasets/wng_main_clean.txt", header = T)

#remove insecticide
main <- main[main$CODE != "WG1P6", ]

# subset tree dataset
tree <- main %>%
  filter(LIFE.FORM %in% c("tree", "shrub"))

# Descriptive statistics
test <- read.table("datasets/wng_main_test.txt")
tree_test <- read.table("datasets/wng_tree_test.txt")

# 0.b Remove insecticide ----
test <- test[rownames(test) != "WG1P6", ]
tree_test <- tree_test[rownames(tree_test) != "WG1P6", ]
ggplot(test, aes(y=log(BIO), x=TREAT)) + geom_point()


# Statistical models -----

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

# >>> FIG. 1 -----
fig1 <- ggplot(panel_plot_staked_data, aes(x=treat, y=val)) + 
  facet_wrap(~type, scales = "free") + 
  geom_point(col = "grey80", alpha = 0.5, cex = 4) + 
  theme_bw()

colors = color=c("black","grey50","grey50","grey50","grey50","red",
                 "black","grey50","red","grey50","grey50","grey50",
                 "black","red","grey50","grey50","grey50","red")

#png("figs/fig1.png",width=1200, height = 400)
fig1 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color= colors, width=0.2, lwd=1.5) +
  stat_summary(fun.y=mean, geom="point", color=colors, cex = 5) +
  theme(axis.text.x=element_text(angle=0, size=20, hjust=0.5),
        axis.text.y=element_text(angle=0, size=20, hjust=0.5),
      strip.text = element_text(size=20),
      legend.justification=c(0.5,0.5), 
      legend.position="bottom")+xlab("")+ylab("")
#dev.off()

# >>> FIG. 2 ----

# cwm_analysis.R

# >>> FIG. 3 SPECIES LEVEL ANALYSIS----

# make subsets

# main and tree datasets
comb <- expand.grid(unique(main$TREAT)[3], unique(main$TREAT))
comb <- comb[-which(comb[,1] == comb[,2]),]

# I will focus only on the tree species as suggsted by ordination graph

# for (comb_no in 1:5){
#  
# }


# All data sub_set!
cvst <- list() # for most abundant species comparisons
unique_spec <- list() # unique species datasets
pair_comp <- list() # community comparisons RDA - trees
pair_comp_treat <- 

for (comb_no in 1:dim(comb)[1]){
  # step 1 subset the dataset
  subset <- tree[tree$TREAT %in% c(as.character(comb[comb_no,]$Var1), 
                                   as.character(comb[comb_no,]$Var2)), ]
  
  # subset_main <- main[main$TREAT %in% c(as.character(comb[comb_no,]$Var1), 
  #                                  as.character(comb[comb_no,]$Var2)), ]
  
  usp_name_list <- paste(comb[comb_no,]$Var1,
                         comb[comb_no,]$Var2, sep="_")
  
  pair_comp[[usp_name_list]] <- contingencyTable2(subset, "SP_CODE", "CODE","WEIGHT")
  
  # subset_RDA_main <- contingencyTable2(subset_main, "SP_CODE", "CODE","WEIGHT")
  
  # step 2 most often occuring species present in both
  sub_table <- table(subset$SP_CODE, subset$TREAT)
  
  # Species present in a control plot but not in the treatment
  sp_present <- sub_table[,c(comb[comb_no,]$Var1,comb[comb_no,]$Var2)]
  sp_present <- sp_present[rowSums(sp_present) != 0, ]
  
  entry <- list(unique_treat = which(sp_present[,1] == 0),# zero at treatment
       unique_con = which(sp_present[,2] == 0)) # zero at control
  
  unique_spec[[usp_name_list]] <- entry
  #tree_ct <- contingencyTable2(tree, "SP_CODE", "TREAT", "WEIGHT")
  
  # Which species occur at leas "threshold" times
  threshold <- 7
  selected <- names(which(rowSums(sub_table) >= threshold))
  selected_subset <- subset[subset$SP_CODE %in% selected, ]
  
  # Dataset list C vs I
  for (spec in selected){
    # Plot for a given species: Piptar
    sub_dat <- selected_subset[selected_subset$SP_CODE == spec, ]
    name_list <- paste(comb[comb_no,]$Var1,
                       comb[comb_no,]$Var2, 
                       spec,sep="_")
    cvsi[[name_list]] <- contingencyTable2(sub_dat, "BLOCK", "TREAT", "WEIGHT")
  }
}


# see if they are most abundant as well
biomass <- tapply(subset$WEIGHT, subset$SP_CODE,sum, na.rm = TRUE)
biomass[order(biomass, decreasing = T)]

colors <- rep("black", length(selected)*2)

# run tests
# spec <- selected[1]
# for (spec in selected){
#   sp_dataset <- subset[subset$SP_CODE == spec,]
#   #model <- lm(log(WEIGHT) ~ TREAT, data = sp_dataset)
#   # model <- brm(WEIGHT ~ TREAT + (1|BLOCK), 
#   #              family=lognormal(),
#   #              data=sp_dataset)
#   # 
#   # fv <- fitted(model, summary = FALSE)
#   # m <- mean(fv)
#   # prob_greater <- colMeans(fv > m)
#   # prob_smaller <- 1 - prob_greater
#   
#   print(summary(model))
#   stanplot(model, type = "hist")
# }

# make plots
panelspecies <- ggplot(selected_subset, aes(x = TREAT, y = WEIGHT)) +
  facet_wrap(~SP_CODE, scales = "free") + 
  geom_point(col = "grey80", alpha = 0.5, cex = 4) + 
  geom_line(aes(linetype = BLOCK), size = 0.5, alpha=0.15) +
  theme_bw()

panelspecies + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                            geom="errorbar", color= colors, width=0.2, lwd=1.5) +
  stat_summary(fun.y=mean, geom="point", color= colors, cex = 5) +
  theme(axis.text.x=element_text(angle=0, size=20, hjust=0.5),
        axis.text.y=element_text(angle=0, size=20, hjust=0.5),
        strip.text = element_text(size=20),
        legend.justification=c(0.5,0.5), 
        legend.position="bottom")+xlab("")+ylab("")

# See the number of observations
selected_subset$SP_CODE <- as.character(selected_subset$SP_CODE)
selected_subset$TREAT <- as.character(selected_subset$TREAT)
table(selected_subset$SP_CODE, selected_subset$TREAT)
        