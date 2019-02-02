# 0. Load data and packages ----

library(dplyr)
library(lme4)
library(lmerTest)
library(Hmisc) #for the stat_summary plots
library(brms)
library(car)
library(brms)
library(vegan)
library(MASS)

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


# Statistical models -----

# Cumulative biomass
bio_rbl <- lmer(log(BIO) ~ TREAT + (1|GARDEN),
                data = test) 

bio_rbl_t <- lmer(logBio ~ TREATMENT + (1|GARDEN),
                data = tree_test) 

summary(bio_rbl)
summary(bio_rbl_t)

# Richness
spn_rbl <- lmer(SPEC_NO ~ TREAT + (1|GARDEN),
                data = test)

spn_rbl2 <- glmer(SPEC_NO ~ TREAT + (1|GARDEN),
                 data = test, family = "poisson")


AIC(spn_rbl, spn_rbl2)
summary(spn_rbl2)
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

# main and tree datasets
comb <- expand.grid(unique(main$TREAT)[3], unique(main$TREAT))
comb <- comb[-which(comb[,1] == comb[,2]),]

# >> data summary -----
# I will focus only on the tree species as suggsted by ordination graph
# All data subsets. Analyse relative proportion of each species in the community,
# and not the absolute biomass. In that case all of these will simply decrease.

cvst <- list()            # for most abundant species comparisons
unique_spec <- list()     # unique species datasets
pair_comp <- list()       # community comparisons RDA - trees
pair_comp_treat <- list() # treatments

# Relative abundance data
ratree <- tree[order(tree$CODE), ]
ratree <- main[order(main$CODE), ] # for the whole community
ratree$WEIGHT <-  stack(tapply(tree$WEIGHT, tree$CODE, function(x){x/sum(x)}))$value
ratree$WEIGHT <-  stack(tapply(main$WEIGHT, main$CODE, function(x){x/sum(x)}))$value

# test (any plot code should show same species)
# code <- "WG2P2"
# a <- ratree[ratree$CODE == code, c("WEIGHT","SP_CODE")]
# b <-   tree[  tree$CODE == code, c("WEIGHT","SP_CODE")]
# plot(a$WEIGHT~b$WEIGHT)

#comb_no <- 2
for (comb_no in 1:dim(comb)[1]){
  # step 1 subset the dataset to a given pair of treatments
  subset <- ratree[ratree$TREAT %in% c(as.character(comb[comb_no,]$Var1), #control
                                   as.character(comb[comb_no,]$Var2)), ]
  
  # subset_main <- main[main$TREAT %in% c(as.character(comb[comb_no,]$Var1), 
  #                                  as.character(comb[comb_no,]$Var2)), ]
  
  usp_name_list <- paste(comb[comb_no,]$Var1,
                         comb[comb_no,]$Var2, sep="_")
  
  ct_sub <- contingencyTable2(subset, "CODE","SP_CODE", "WEIGHT")
  rowSums(ct_sub)
  
  pair_comp[[usp_name_list]] <- (ct_sub)
  
  treats <- subset[, c("CODE", "BLOCK","TREAT")]
  pair_comp_treat[[usp_name_list]] <- treats[!duplicated(treats),] #All plots with their treatments.
  
  # convert weights into relative weights
  
  # step 2 most often occuring species present in both
  sub_table <- table(subset$SP_CODE, subset$TREAT)
  
  # Species present in a control plot but not in the treatment
  sp_present <- sub_table[,c(comb[comb_no,]$Var1,comb[comb_no,]$Var2)]
  sp_present <- sp_present[rowSums(sp_present) != 0, ]
  
  entry <- list(which(sp_present[,1] == 0), # zero at treatment
                which(sp_present[,2] == 0)) # zero at control
  
  unique_species_biomass1 <- pair_comp[[usp_name_list]][, names(entry[[1]])]
  unique_species_biomass2 <- pair_comp[[usp_name_list]][, names(entry[[2]])]
  
  u_entry <- list(unique_species_biomass1,
                  unique_species_biomass2)
  
  names(u_entry) <- c(colnames(sp_present)[1],colnames(sp_present)[2])
  
  unique_spec[[usp_name_list]] <- u_entry
  
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
    cvst[[name_list]] <- contingencyTable2(sub_dat, "BLOCK", "TREAT", "WEIGHT")
  }
}

# Remove tremor from cvst data
cvst <- cvst[names(cvst) != "CONTROL_FUNGICIDE_TREMOR"]
length(cvst)
# models and plots for species

stacked_cvst <- data.frame()

for (i in 1:length(names(cvst))){
  nm <- names(cvst)[i]
  splitted <- strsplit(nm, "_") 
  str1 <- substr(splitted[[1]][1], 1, 1)
  str2 <- substr(splitted[[1]][2], 1, 1)
  if (splitted[[1]][2] == "WEEVIL125"){str2 <- "H2"}
  if (splitted[[1]][2] == "WEEVIL25"){str2 <- "H1"}
  comb <- paste(str1, str2,
                sep=" vs ")
  spec <- splitted[[1]][3]
  sp_data <- cvst[[nm]]
  #sp_data <- 1000*cvst[[nm]]
  sp_stack <- stack(as.data.frame(sp_data))
  sp_stack$block <- rownames(sp_data)
  sp_stack$type <- nm
  sp_stack$comb <- comb
  sp_stack$spec <- spec
  stacked_cvst <- rbind(stacked_cvst, sp_stack)
}

stacked_cvst

#get rid of TREMOR form the comparisons
stacked_cvst <- stacked_cvst[stacked_cvst$type != "CONTROL_FUNGICIDE_TREMOR", ]

# Statistical tests for biomass increase for individual species

# Tweedie distriobution
# library(cplm)
# colours <- c()

# # for (name in 1:length(names(cvst))){
#   sp_data <- cvst[[name]]
#   
#   # Add one and se
#   
#   print(names(cvst)[[name]])
#   #sp_data <- 1000*cvst[[nm]]
#   sp_stack <- stack(as.data.frame(sp_data))
#   sp_stack$block <- rownames(sp_data)
# 
# f0 <- cpglmm(log(values+1)~(1|block),
#              sigma = ind, data = sp_stack)
#   f1 <- cpglmm(log(values+1)~ind+ (1|block),
#                data = sp_stack)
#   pval <- anova(f0, f1)$`Pr(>Chisq)`[2]
#   print(pval)
#   colours <- c(colours, "black")
#   if (pval <= 0.05){
#     colours <- c(colours, "red")
#   } else {
#     colours <- c(colours, "grey30")
#       }
# # }


# using classic mixed effect models
colours <- c()

# Order cvst alphabetically
cvst <- cvst[order(names(cvst))]

# par(mfrow = c(1,1))

diff_dat <- data.frame()
  
# name <- 19
for (name in 1:length(names(cvst))){
  sp_data <- cvst[[name]]
  print(names(cvst)[[name]])
  #sp_data <- 1000*cvst[[nm]]
  sp_stack <- stack(as.data.frame(sp_data))
  sp_stack$block <- rownames(sp_data)
  
  #treatment - control
  dd <- data.frame(diftc = c(sp_data[,2] - sp_data[,1]), comp = names(cvst)[[name]])
  diff_dat <- rbind(diff_dat, dd)
  
  # dd <- data.frame( nm = c(sp_data[,2] - sp_data[,1]))
  # colnames(dd) <- names(cvst)[[name]]
  # diff_dat <- cbind(diff_dat, dd)

  
  # If I want to use zero inflated beta, should I also 
  # add blocks where species wasn't present? If it is not present in 
  # a block then there is no information from that? Or is it?
  
  # m0 <- brm(values~(1|block), family=zero_inflated_beta(),
  #            data = sp_stack)
  
  # Zero inflated beta
  # m1 <- brm(values~ ind + (1|block), family=zero_inflated_beta(),
  #           data = sp_stack)
  
  # lt <- sp_stack
  # lt$values <- log((lt$values)/(1-lt$values))
  # lt$values
  # 
  # compare_ic(waic(m0), waic(m1), ic = "waic")
  
  # pred1 <- sp_stack %>% 
  #   tidybayes::add_predicted_draws(m1, n = 1000) %>%
  #   filter(ind == unique(sp_stack$ind)[1])
  # 
  # pred2 <- sp_stack %>% 
  #   tidybayes::add_predicted_draws(m1, n = 1000) %>%
  #   filter(ind == unique(sp_stack$ind)[2])
  
  # Difference between posterior distributions
  
  # hist(pred1$.prediction, breaks = 150)
  # hist(pred2$.prediction, breaks = 150, add=T)
  # 
  # sum(pred1$.prediction - pred2$.prediction < 0)/length(pred1$.prediction)
  # 
  # sum((pred1$.prediction - pred2$.prediction) > 0.05)/length(pred1$.prediction)
  
  # sph <- stanplot(m1, type = "hist", pars="^b_")
  # sph$data
  # 
  # sp_stack %>% 
  #   tidybayes::add_predicted_draws(m1, n = 1000) %>%
  #   ggplot(aes(y = logit(.prediction), x = ind, group = block, fill=block, color = block)) + 
  #   geom_jitter(alpha = 0.5) + 
  #   geom_point(aes(y = logit(values), x = ind), size = 3, alpha = 0.1,
  #              color = "grey50") + 
  #   theme_bw()
  
  # Test for heteroscedasticity
  
  # Differences and paired-t-test
  
  # ltest <- leveneTest(values~ind, data=sp_stack)
  # 
  # print(ltest)
  # 
  # qqnorm(sp_stack$values)
  # qqline(sp_stack$values)
  # 
  # Paired t test
  res <- t.test(sp_data[,1],sp_data[,2],paired = TRUE, alternative = "two.sided")
  pval <- res$p.value
  print(pval)
  
  # Sined rank test
  # res <- wilcox.test(sp_data[,1],sp_data[,2], paired=TRUE)
  # pval <- res$p.value 
  # print(pval)
  
  # Paired t test logit
  # sp_logit <- log(sp_data/(1-sp_data))
  # sp_logit[sp_logit == -Inf] <- 0
  # sp_logit[sp_logit == Inf] <- 0
  # 
  # #boxplot(sp_logit, main=names(cvst)[[name]])
  # 
  # res <- t.test(sp_logit[,1],sp_logit[,2],paired = TRUE, alternative = "two.sided")
  # 
  # pval <- res$p.value
  # print(pval)
  
  # qqnorm(sp_logit)
  # qqline(sp_logit)
  # 
  # valuesBCMod <- caret::BoxCoxTrans(sp_stack$values+1)
  # sp_stack$values <- predict(valuesBCMod, 
  #                            sp_stack$values)
  
  # f0 <- lmer(log(values+1)~(1|block),
  #            data = sp_stack)
  # f1 <- lmer(log(values+1)~ind+ (1|block),
  #              data = sp_stack)
  # pval <- anova(f0, f1)$`Pr(>Chisq)`[2]
  # print(pval)
  colours <- c(colours, "black")
  if (pval <= 0.1){
    colours <- c(colours, "orange")
  }
  if (pval <= 0.05){
    colours <- c(colours, "red")
  } else {
    colours <- c(colours, "grey30")
  }
}

# difference correlations
# diff_dat
# 
# un <- unique(diff_dat$comp)
# d_comb_i <- expand.grid(un[grep("INSECTICIDE", un)], un[grep("INSECTICIDE", un)])
# d_comb_i <- d_comb_i[d_comb_i$Var1 != d_comb_i$Var2, ]
# 
# i <- 1
# 
# for(i in 1:nrow(d_comb_i)){
#   one <- diff_dat[diff_dat$comp == d_comb_i[i,1], ]
#   two <- diff_dat[diff_dat$comp == d_comb_i[i,2], ]
# }

# >>> FIG. 3. B using differences ----

# under construction
# colours <- c()
# 
# for (name in 1:length(names(cvst))){
#   sp_data <- cvst[[name]]
#   print(names(cvst)[[name]])
#   
#   diffs <- sp_data[,1]-sp_data[,2]
#   
#   print(t.test(diffs)$p.value)
  # Stack up the data
  # #sp_data <- 1000*cvst[[nm]]
  # sp_stack <- stack(as.data.frame(sp_data))
  # sp_stack$block <- rownames(sp_data)
  # 
  # log(sp_stack[sp_stack$ind == ,]$values + 1) - log(sp_stack[6:10,]$values + 1)
  # 
  # f0 <- lmer(log(values+1)~(1|block),
  #            data = sp_stack)
  # f1 <- lmer(log(values+1)~ind+ (1|block),
  #            data = sp_stack)
  # pval <- anova(f0, f1)$`Pr(>Chisq)`[2]
  # print(pval)
  # colours <- c(colours, "black")
  # if (pval <= 0.05){
  #   colours <- c(colours, "red")
  # } else {
  #   colours <- c(colours, "grey30")
  # }
# }

# actuall plot ------

# lines log +1
# p1 <- ggplot(stacked_cvst, aes(x = ind, y= log(values+1), group=block)) +
#   geom_line(aes(linetype = block), size = 0.5, alpha=0.15) +
#   facet_wrap(~type, scales = "free", nrow=5, ncol =4) +
#   geom_point(size = 1.7) + theme_bw()

# raw

# p1 <- ggplot(stacked_cvst, aes(x = ind, y= values, group=block)) +
#   geom_line(aes(linetype = block), size = 1, alpha=0.15) +
#   facet_wrap(comb ~ spec, scales = "free", nrow=5, ncol =4) +
#   geom_point(cex = 4) + theme_bw() + 
#   theme(axis.text.x=element_text(angle=0, size=20, hjust=0.5),
#         axis.text.y=element_text(angle=0, size=20, hjust=0.5),
#         strip.text = element_text(size=20),
#         legend.justification=c(0.5,0.5), 
#         legend.position="bottom")+
#   xlab("") + 
#   ylab("")
# p1

p1 <- ggplot(stacked_cvst, aes(x = ind, y= values, group=block)) +
  geom_line(aes(linetype = block), size = 1, alpha=0.15) +
  facet_grid(spec ~ comb, scales = "free") +
  geom_point(cex = 2) + theme_bw() + 
  theme(axis.text.x=element_text(angle=0, size=10, hjust=0.5),
        axis.text.y=element_text(angle=0, size=10, hjust=0.5),
        strip.text = element_text(size=10),
        legend.justification=c(0.5,0.5), 
        legend.position="bottom")+
  xlab("") + 
  ylab("")
p1

# Siimple lines
# png("figs/fig3b.png",width=1200, height = 1200)
p1 <- ggplot(stacked_cvst, aes(x = ind, y= values, group=block)) +
  geom_line(size = 2, alpha=0.15) +
  facet_grid(spec ~ comb, scales = "free") +
  geom_point(cex = 3) + theme_bw() + 
  theme(axis.text.x=element_text(angle=0, size=10, hjust=0.5),
        axis.text.y=element_text(angle=0, size=10, hjust=0.5),
        strip.text = element_text(size=20),
        legend.justification=c(0.5,0.5), 
        legend.position="bottom")+
  xlab("") + 
  ylab("")
p1
# dev.off()

#png("figs/fig3.png",width=1200, height = 1200)
p1 <- ggplot(stacked_cvst, aes(x = ind, y= values, group=block)) +
  geom_line(aes(linetype = block), size = 3, alpha=0.15) +
  facet_grid(spec ~ comb, scales = "free") +
  geom_point(cex = 6) + theme_bw() + 
  theme(axis.text.x=element_text(angle=0, size=20, hjust=0.5),
        axis.text.y=element_text(angle=0, size=20, hjust=0.5),
        strip.text = element_text(size=20),
        legend.justification=c(0.5,0.5), 
        legend.position="bottom")+
  xlab("") + 
  ylab("")
p1
#dev.off()

# summaries, but these go below zero!

# png("figs/fig3sig.png",width=800, height = 800)
p1 <- ggplot(stacked_cvst, aes(x = ind, y= values)) +
  facet_wrap(~type, scales = "free", ncol=4, nrow=5) +
  geom_point(size = 3, col = "grey80") + theme_bw()

p1 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                geom="errorbar", color =colours, width=0.2, lwd=1.5) +
  stat_summary(fun.y=mean, geom="point", color =colours, cex = 5)
# dev.off()

stacked_cvst_logit <- stacked_cvst
stacked_cvst_logit$values <- log(stacked_cvst_logit$values/(1-stacked_cvst_logit$values))
stacked_cvst_logit$values[stacked_cvst_logit$values == -Inf] <- 0

qqnorm(stacked_cvst_logit$values)
qqline(stacked_cvst_logit$values)

p1 <- ggplot(stacked_cvst_logit, aes(x = ind, y= values)) +
  facet_wrap(~type, scales = "free", ncol=4, nrow=5) +
  geom_point(size = 3, col = "grey80") + theme_bw()

p1 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                  geom="errorbar", color =colours, width=0.2, lwd=1.5) +
  stat_summary(fun.y=mean, geom="point", color =colours, cex = 5)


# some good models for this?
# tweedie continuous data?

# example! 

qqnorm(log(sp_stack[1:5,]$values + 1) - log(sp_stack[6:10,]$values + 1))
qqline(log(sp_stack[1:5,]$values + 1) - log(sp_stack[6:10,]$values + 1))
(sp_stack[1:5,]$values)
qqnorm(sp_stack[6:10,]$values)



sp_stack
# bayesian
f2 <- bcplm(values~ind+(1|block),
            data = sp_stack, n.iter = 15000,
              n.burnin = 5000, n.thin = 10)

# gelman.diag(f2$sims.list)
# summary(f2)
# 
# acfplot(f2$sims.list, lag.max = 20)
# xyplot(f2$sims.list)                              
# densityplot(f2$sims.list)               
# summary(f2)
# plot(f2)
# 
# library(nlme)
# l1 <- lme(log(values+1)~ind, random = ~1|block,
#           data = sp_stack)
# l2 <- lme(log(values+1)~ind, random = ~1|block,
#     weights=varIdent(form=~1|ind),
#     data = sp_stack)
# 
# AIC(l1, l2)
# 
# summary(l1)
# nod1 <- brm(log(values+1)~ind+(1|block), data = sp_stack)
# 
# # plots for the bayesian model
# sp_stack %>% 
#   tidybayes::add_predicted_draws(nod1, n = 1000) %>%
#   ggplot(aes(y = log(.prediction+1), x = ind)) + 
#   geom_jitter() +
#   theme_bw()



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

## >>> FIG X. RDA once again -----

# instead of kg If I use grams, doesnt matter. What I want are the proportions!!!
# relative distributions!

# Rda with relative abundances
ct_main <- contingencyTable2(main, "CODE","SP_CODE", "WEIGHT")
ra <- ct_main/rowSums(ct_main)

# make dataset only for trees
ct_tree <- contingencyTable2(tree, "CODE","SP_CODE", "WEIGHT")
ra_tree <- ct_tree/rowSums(ct_tree)

# Order and rename the test dataset
test <- test[rownames(ra), ]
names(test) <- c("SPEC_NO","SW","SIMP","EVEN","BIO","GARDEN","PLOT","PLOT_CODE","TREAT",
                 "LOCATION","control","fungicide","herbivory_moderate","herbivory_high",
                 "insecticide","predator_exclosure","DIST")

# RDA models -----
m <- rda(ra~fungicide+herbivory_moderate+herbivory_high+insecticide+predator_exclosure+Condition(GARDEN),
         data=test, scale=FALSE)

# m_ca <- cca(ra)
# plot(m_ca, display = "sites")
# m <- rda(relative_abundances~FUN+HLO+HHI+INS+PRE+Condition(GARDEN),data=test)
anova(m, by="terms", permutations = 999)
anova(m, by="axis", permutations = 999)

summary(m)

coef(m)
R2 <- RsquareAdj(m)$r.squared
R2 
 
# Check the effects for trees
test <- test[rownames(ra_tree), ]
mt <- rda(ra_tree~fungicide+herbivory_moderate+herbivory_high+insecticide+predator_exclosure+Condition(GARDEN),
         data=test)
# m <- rda(relative_abundances~FUN+HLO+HHI+INS+PRE+Condition(GARDEN),data=test)
anova(mt, by="terms", permutations = 9999)
# plot(mt)

# # goodness of fit ----
gof <- goodness(m, display = "species", statistic = "explained",
         summarize = FALSE)
gof2 <- inertcomp(m, display = "species", proportional = TRUE)

specs <-sort(round(gof2[,1], 3), decreasing = T)
selected <- names(specs[1:10])

specs <-sort(round(gof[,1], 3), decreasing = T)
selected <- names(specs[1:10])



# actuall figure ----

#https://www.fromthebottomoftheheap.net/2012/04/11/customising-vegans-ordination-plots/

scl = 2
#plot(m, scaling = scl)
plot(m, type = "n", scaling = scl)
colvec <- c(rgb(230,159,0, max=255), 
            rgb(86,180,233, max=255), 
            rgb(0,158,115, max=255), 
            rgb(0,114,178, max=255), 
            rgb(213,94,0, max=255), 
            rgb(240,228,66, max=255))
with(test, points(m, display = "sites", col = colvec[TREAT],
                      scaling = scl, pch = 21, cex = 2, bg = colvec[TREAT]))
# Species, but it is just a mess
# text(m, display = "species", scaling = scl, cex = 0.8,
#      col = rgb(1,1,1, max = 255, alpha=100))

# Legend
with(test, legend("topright", legend = levels(TREAT), bty = "n",
                      col = colvec, pch = 21, pt.bg = colvec))
text(m, display="bp", col="darkcyan", lwd=2)

selection <- rownames(summary(m)$species) %in% selected
orditorp(m, scaling = -1, 
         display = "species", select = selection, col = "darkcyan")

# orditkplot(m,scaling=scl, display = "species")

# variance explained ----
barplot(m$CA$eig/m$tot.chi, names.arg = 1:m$CA$rank, cex.names = 0.5, 
        ylab="Proportion of variance explained", xlab="CA axis")
m$CCA$tot.chi/m$tot.chi



# >> Pairwise comparisons of the community structure RDA -----
ds_num <- 5

dts <- t(pair_comp[[ds_num]])
radts <- dts/rowSums(dts)

# Order the test dataset
ratest <- pair_comp_treat[[ds_num]]
rownames(ratest) = ratest$CODE
ratest <- ratest[rownames(radts), ]

m1 <- rda(log(radts+1)~TREAT+Condition(BLOCK), data = ratest)
anova(m1, by="terms", permutations = 9999)
plot(m1)

# >> dominance plots ----
radts <- log(radts+1) # logged values
# ds <- log(1000*pair_comp[[ds_num]]+1) # logged grams
# ds <- pair_comp[[ds_num]] # raw

trt_a <-unique(ratest$TREAT)[1]
trt_b <-unique(ratest$TREAT)[2]

ds_a <- radts[ratest$TREAT == trt_a, ]
ds_b <- radts[ratest$TREAT == trt_b, ]

par(mfrow = c(1,2))

ds_a[ds_a == 0] <- NA
lineup <- order(colSums(ds_a, na.rm = T))
boxplot(ds_a[, rev(lineup)], main = trt_a, las=2)

ds_b[ds_b == 0] <- NA
lineup <- order(colSums(ds_b, na.rm = T))
boxplot(ds_b[, rev(lineup)], main = trt_b, las=2)






### experiments with hurdle lognormal model

# vignette:
#https://cloud.r-project.org/web/packages/GLMMadaptive/vignettes/ZeroInflated_and_TwoPart_Models.html

# library(GLMMadaptive)
# 
# 
# set.seed(1234)
# n <- 100 # number of subjects
# K <- 8 # number of measurements per subject
# t_max <- 5 # maximum follow-up time
# 
# # we constuct a data frame with the design: 
# # everyone has a baseline measurment, and then measurements at random follow-up times
# DF <- data.frame(id = rep(seq_len(n), each = K),
#                  time = c(replicate(n, c(0, sort(runif(K - 1, 0, t_max))))),
#                  sex = rep(gl(2, n/2, labels = c("male", "female")), each = K))
# 
# # design matrices for the fixed and random effects non-zero part
# X <- model.matrix(~ sex * time, data = DF)
# Z <- model.matrix(~ time, data = DF)
# # design matrices for the fixed and random effects zero part
# X_zi <- model.matrix(~ sex, data = DF)
# Z_zi <- model.matrix(~ 1, data = DF)
# 
# betas <- c(-2.13, -0.25, 0.24, -0.05) # fixed effects coefficients non-zero part
# sigma <- 0.5 # standard deviation error terms non-zero part
# gammas <- c(-1.5, 0.5) # fixed effects coefficients zero part
# D11 <- 0.5 # variance of random intercepts non-zero part
# D22 <- 0.1 # variance of random slopes non-zero part
# D33 <- 0.4 # variance of random intercepts zero part
# 
# # we simulate random effects
# b <- cbind(rnorm(n, sd = sqrt(D11)), rnorm(n, sd = sqrt(D22)), rnorm(n, sd = sqrt(D33)))
# # linear predictor non-zero part
# eta_y <- as.vector(X %*% betas + rowSums(Z * b[DF$id, 1:2, drop = FALSE]))
# # linear predictor zero part
# eta_zi <- as.vector(X_zi %*% gammas + rowSums(Z_zi * b[DF$id, 3, drop = FALSE]))
# # we simulate log-normal longitudinal data
# DF$y <- exp(rnorm(n * K, mean = eta_y, sd = sigma))
# # we set the zeros from the logistic regression
# DF$y[as.logical(rbinom(n * K, size = 1, prob = plogis(eta_zi)))] <- 0
# 
# hurdle.lognormal <- function () {
#   stats <- make.link("identity")
#   log_dens <- function (y, eta, mu_fun, phis, eta_zi) {
#     sigma <- exp(phis)
#     # binary indicator for y > 0
#     ind <- y > 0
#     # non-zero part
#     eta <- as.matrix(eta)
#     eta_zi <- as.matrix(eta_zi)
#     out <- eta
#     out[ind, ] <- plogis(eta_zi[ind, ], lower.tail = FALSE, log.p = TRUE) + 
#       dnorm(x = log(y[ind]), mean = eta[ind, ], sd = sigma, log = TRUE)
#     # zero part
#     out[!ind, ] <- plogis(eta_zi[!ind, ], log.p = TRUE)
#     attr(out, "mu_y") <- eta
#     out
#   }
#   score_eta_fun <- function (y, mu, phis, eta_zi) {
#     sigma <- exp(phis)
#     # binary indicator for y > 0
#     ind <- y > 0
#     # non-zero part
#     eta <- as.matrix(mu)
#     out <- eta
#     out[!ind, ] <- 0
#     out[ind, ] <- (log(y[ind]) - eta[ind, ]) / sigma^2
#     out
#   }
#   score_eta_zi_fun <- function (y, mu, phis, eta_zi) {
#     ind <- y > 0
#     probs <- plogis(as.matrix(eta_zi))
#     out <- 1 - probs
#     out[ind, ] <- - probs[ind, ]
#     out
#   }
#   score_phis_fun <- function (y, mu, phis, eta_zi) {
#     sigma <- exp(phis)
#     # binary indicator for y > 0
#     ind <- y > 0
#     # non-zero part
#     eta <- as.matrix(mu)
#     out <- eta
#     out[!ind, ] <- 0
#     out[ind, ] <- - 1 + (log(y[ind]) - eta[ind, ])^2 / sigma^2
#     out
#   }
#   simulate <- function (n, mu, phis, eta_zi) {
#     y <- rnorm(n = n, mean = mu, sd = exp(phis))
#     y[as.logical(rbinom(n, 1, plogis(eta_zi)))] <- 0
#     y
#   }
#   structure(list(family = "two-part log-normal", link = stats$name, 
#                  linkfun = stats$linkfun, linkinv = stats$linkinv, log_dens = log_dens,
#                  score_eta_fun = score_eta_fun, score_eta_zi_fun = score_eta_zi_fun,
#                  score_phis_fun = score_phis_fun, simulate = simulate),
#             class = "family")
# }
# 
# km1 <- mixed_model(y ~ sex * time, random = ~ 1 | id, data = DF, 
#                    family = hurdle.lognormal(), n_phis = 1,
#                    zi_fixed = ~ sex)
# 
# km1
# 
# km1 <- mixed_model(values ~ ind, random = ~ 1 | block, data = sp_stack, 
#                    family = hurdle.lognormal(), n_phis = 1)


#https://stats.stackexchange.com/questions/339610/hurdle-model-with-non-zero-gaussian-distribution-in-r
#https://stats.idre.ucla.edu/r/dae/tobit-models/

# Description of the comunity ----

# most abundant species
compare <- expand.grid(as.character(unique(main$TREAT)),
                       as.character(unique(main$TREAT)))
compare <- compare[compare$Var1 != compare$Var2, ]
compare <- compare[compare$Var1 == "CONTROL", ]
compare[,1] <- as.character(compare[,1])
compare[,2] <- as.character(compare[,2])


mainC <- main[main$TREAT == "CONTROL", ]
library(forcats)


# strings <- c("cf", "ch2", "cp", "ch1", "ci")
# 
# for (i in 1:dim(compare)[1]){
#   cnms <- c(compare[i,]$Var1, compare[i,]$Var2)
#   dst <- main[main$TREAT %in% c(compare[i,]$Var1, compare[i,]$Var2), ]
#   path <- paste("figs/figS", strings[i] ,".png", sep = "")
#   png(path,width=1200, height = 600)
#   p <- ggplot(dst, aes(x = fct_reorder(SP_CODE, WEIGHT, fun = median, .desc =TRUE), 
#                        y = log(WEIGHT), fill=TREAT)) + 
#     geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 90, 
#                                                                    hjust = 1,
#                                                                    vjust = 0.5,
#                                                                    size = 10)) +
#     scale_x_discrete(name="") +
#     scale_y_continuous(name="log[kg]") + 
#     ggtitle("Dominance structure")
#   
#   print(p)
#   dev.off()
#   readline(prompt="Press [enter] to continue")
# }

# Melochia sp1 biomass
exp(4.85957982)

# Average biomass per plot:
w <- as.vector(tapply(mainC$WEIGHT, mainC$CODE, sum, na.rm = TRUE))
w <- w[!is.na(w)]
summary(w)

library(EnvStats)
vals <- elnorm(w, ci=TRUE)
names(vals)
exp(vals$parameters[1])
exp(3.846207)
exp(4.971028)

# Average BC values
sp_data_c <- contingencyTable2(mainC, "SP_CODE", "CODE", "WEIGHT")
summary(vegdist(t(sp_data_c), method="bray"))

treeC <- tree[tree$TREAT == "CONTROL", ]
sp_tree_c <- contingencyTable2(treeC, "SP_CODE", "CODE", "WEIGHT")
mean(vegdist(t(sp_tree_c), method="bray"))

# TRAITS for most abundant species ----

#Piptar," mikamia", "melamu"
pmds <- main[main$SP_CODE %in% c("PIPTAR", "MIKAMI", "MELAMU", 
                                 "PIPEUM", "TRICPL"), ]
pmds <- pmds[pmds$TREAT == "CONTROL", ]
pmds$SP_CODE <- as.character(pmds$SP_CODE) 
boxplot(HERB~SP_CODE, data=pmds)

# inverse logistic
il <- function(x){
  return(exp(x)/(exp(x)+1))
}

library(emmeans)

# herbBeta <- brm(il(HERB)~SP_CODE, data = pmds, family = Beta())
# stanplot(herbBeta, pars="^b_")
ggplot(pmds, aes(y=il(HERB), x=SP_CODE)) + geom_point()
#leveneTest(HERB~SP_CODE, data=pmds, family=Beta()) # heteroscedasticity
# glmmTMB(il(HERB) ~ SP_CODE + (1|BLOCK), data=pmds, family=list(family="beta",link="logit")) 
lmerherb <- lmer(HERB~SP_CODE+(1|BLOCK), data=pmds)
phtherb <- emmeans(lmerherb, list(pairwise~SP_CODE), adjust = "tukey")
#plot(phtherb)
p1 <- emmip(lmerherb, ~ SP_CODE, CIs = TRUE) + ggtitle("Herbivory")+ 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

ggplot(pmds, aes(y=SLA, x=SP_CODE)) + geom_point()
#leveneTest(SLA~SP_CODE, data=pmds)  # no heteroscedasticity
lmersla <- lmer(SLA~SP_CODE+(1|BLOCK), data=pmds)
phtsla <- emmeans(lmersla, list(pairwise~SP_CODE), adjust = "tukey")
#plot(phtsla)
p2 <- emmip(lmersla, ~ SP_CODE, CIs = TRUE) + ggtitle("SLA")+ 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

# LDMC tukey
ggplot(pmds, aes(y=LDMC, x=SP_CODE)) + geom_point()
#leveneTest(LDMC~SP_CODE, data=pmds)  # no heteroscedasticity
lmerldmc <- lmer(LDMC~SP_CODE+(1|BLOCK), data=pmds)
phtldmc <- emmeans(lmerldmc, list(pairwise~SP_CODE), adjust = "tukey")
#plot(phtldmc)
p3 <- emmip(lmerldmc, ~ SP_CODE, CIs = TRUE) + ggtitle("LDMC") + 
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# png("figs/figXtukey.png", width=800, height = 300)
multiplot(p1, p2, p3, cols=3)
# dev.off()

# Trying with ordisurf
source("code/cwm_analysis.R")
traitdat <- traitdat[rownames(m$CA$u), ]
par(mfrow=c(1,3))
sherb <- ordisurf(m, traitdat$cwmherb,family="gaussian", col="darkcyan",
                  main = "Herbivory")
ssla  <- ordisurf(m, traitdat$cwmsla,family="gaussian", col="darkcyan",
                  main = "SLA")
sldmc <- ordisurf(m, traitdat$cwmldmc,family="gaussian", col="darkcyan",
                  main = "LDMC")

# All significant??
par(mfrow=c(1,1))
summary(sherb)
summary(ssla)
summary(sldmc)
