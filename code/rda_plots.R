library(ggplot2)
library(vegan)
library(ggord)

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



#write.table(WNGmain, "wng_main.txt")
wng_biomass <- main
WNGtreat <- test

#Drop the 
#Swith off the Insecticide plot WG1P6 from the biomass data
# wng_biomass <- wng_biomass[wng_biomass$CODE != "WG1P6", ]
# WNGtreat <- WNGtreat[WNGtreat$PLOT_CODE != "WG1P6",]

wng_biomass_trees <- wng_biomass[wng_biomass$LIFE.FORM %in% c("tree","shrub"),]
mat_aw <-contingencyTable2(wng_biomass, "CODE","SP_CODE","WEIGHT")
mat_atw <-contingencyTable2(wng_biomass_trees, "CODE","SP_CODE","WEIGHT")
mat_nmds <- mat_aw
# proportions:
mat_aw <- mat_aw/rowSums(mat_aw)

# mat_atw <- read.table("mat_atw.txt")
# proportions
mat_atw <- mat_atw/rowSums(mat_atw)

# Reduced without G1 (ecotone, close to the vilage)
red_mat_aw <- mat_aw[-(1:6),]

WNGtreatRed <- WNGtreat[WNGtreat$GARDEN != "WG1",]

# RDAall <- vegan::rda(log(as.data.frame(mat_aw)+1)~FUN+HLO+HHI+PRE+INS+Condition(GARDEN),
#                      data=WNGtreat[rownames(WNGtreat) %in% rownames(mat_aw),])
# RDAtree <- vegan::rda(log(as.data.frame(mat_atw)+1)~FUN+HLO+HHI+PRE+INS+Condition(GARDEN),
#                       data=WNGtreat[rownames(WNGtreat) %in% rownames(mat_atw),])

# ######################## NMDS of all plots
# data <-contingencyTable2(wng_biomass, "CODE","SP_CODE","WEIGHT")
# WNGmain[WNGmain$SP_CODE == "MIMODI",]
# wng_biomass[wng_biomass$SP_CODE == "MIMODI",]
# 
# vegdist(data)
# 
# allNMDS <- metaMDS(data)
# ordiplot(allNMDS,display="sites")
# orditorp(allNMDS,display="species",col="red",air=0.01)
# orditorp(allNMDS,display="sites",cex=1.25,air=0.01)

RDAall <- vegan::rda(as.data.frame(mat_aw)~FUN+HLO+HHI+PRE+INS+Condition(GARDEN),
                     data=WNGtreat[rownames(WNGtreat) %in% rownames(mat_aw),])
# RDAallred <- vegan::rda(as.data.frame(red_mat_aw)~FUN+HLO+HHI+PRE+INS+Condition(GARDEN),
#                      data=WNGtreatRed[rownames(WNGtreatRed) %in% rownames(red_mat_aw),])

RDAtree <- vegan::rda(mat_atw~FUN+HLO+HHI+PRE+INS+Condition(GARDEN),
                      data=WNGtreat[rownames(WNGtreat) %in% rownames(mat_atw),])

#par(mfrow=c(1,2))

pAll <- ggord(RDAall, as.character(WNGtreat$TREATMENT), ellipse=TRUE,
              ptslab = TRUE,alpha_el = 0.1) + theme(legend.position = 'none')
pAll

# Reduced dataset without G1
# pAllred <- ggord(RDAallred, as.character(WNGtreatRed$TREATMENT), ellipse=TRUE,
#               ptslab = TRUE,alpha_el = 0.1)
# pAllred

pTree <- ggord(RDAtree, as.character(WNGtreat$TREATMENT), ellipse=TRUE,
               ptslab = TRUE,alpha_el = 0.1) + theme(legend.position = 'none')
pTree

# Species responding strongly to the 
# RDA1, related to INS -> FUN
# round(as.matrix(sort(summary(RDAall)$species[,1],decreasing=T)),5)
# # RDA2, No herb -> high HER
# sort(summary(RDAall)$species[,2],decreasing=T)
# 
# round(as.matrix(sort(summary(RDAtree)$species[,1],decreasing=T)),5)
# round(as.matrix(sort(summary(RDAtree)$species[,2],decreasing=T)),5)
# 
# summary(RDAtree)

# Some analyses
#RDAall$CCA$v

#Design based permutations
# h <- how(blocks = WNGtreat$GARDEN, nperm=3000)
# p1 <- anova(RDAall, permutations = h,parallel = 3)
# p1axis <- anova(RDAall, permutations = h,parallel = 3, by="terms")
# p2axis <- anova(RDAtree, permutations = h,parallel = 3, by="terms")
# 
# write.table(p1axis, "rand_res_all_prop.txt")
# write.table(p2axis, "rand_res_tre_prop.txt")
# write.table(p1axis, "rand_res_all_prop_no_1g6.txt")
# write.table(p2axis, "rand_res_tre_prop_no_g1p6.txt")

# setwd("C:/Users/Piotr Szefer/Desktop/Work/garden experiment/datasets")
# # p1axis <- read.table("rand_res_all.txt")
# # p2axis <- read.table("rand_res_tre.txt")
# 
p1axis <- read.table("rand_res_all_prop.txt")
p2axis <- read.table("rand_res_tre_prop.txt")

library(knitr)

######
#NMDS#
######

# nmdstrees <- metaMDS(mat_atw)
# plot(nmdstrees, type ="n")
# colors <- as.numeric(WNGtreat$TREATMENT)
# points(nmdstrees, col=colors, pch=19, cex=1.5)
# text(nmdstrees, display = "sites")
# 
# for(name in as.character(unique(WNGtreat$TREATMENT))){
#   tname <- name
#   rn <- as.character(WNGtreat[WNGtreat$TREATMENT %in% tname,]$PLOT_CODE)
#   vdb <- as.matrix(vegdist(mat_atw))
#   fbc <- vdb[rownames(vdb) %in% rn, colnames(vdb) %in% rn]
#   diag(fbc) <- NA
#   print(paste(name, mean(fbc, na.rm=T)))
# }

################
#Dissimilarity SPECIES#
################

# full_inc <- c()
# for(name in as.character(unique(WNGtreat$TREATMENT))){
#   tname <- name
#   rn <- as.character(WNGtreat[WNGtreat$TREATMENT %in% tname,]$PLOT_CODE)
#   inc_vec <- mat_atw[rownames(mat_atw) %in% rn, ]
#   inc_vec_cs <- colSums(inc_vec)
#   full_inc <- rbind(full_inc, inc_vec_cs)
# }
# rownames(full_inc) <- as.character(unique(WNGtreat$TREATMENT))
# vegdist(full_inc, upper = T)
# 
# 
# full_inc[full_inc > 0] <- 1
# full_inc <- t(full_inc)
# fc <- full_inc[,c(1,3)]
# 
# shared<- length(rownames(fc[(fc[,1] ==1 & fc[,2] == 1),]))
# uniq_f<- length(rownames(fc[(fc[,1] ==1 & fc[,2] != 1),]))
# uniq_c<- length(rownames(fc[(fc[,1] !=1 & fc[,2] == 1),]))

# # 
# # # # Is the position can be explained by traits?
# # 
# sort(inertcomp(RDAtree)[,1], decreasing = T) -> CSTR
# sort(inertcomp(RDAall)[,1], decreasing = T) -> CSTRa
# round(CSTRa, 3)[1:4]
# round(CSTR, 3)[1:4]
# # sort(inertcomp(RDAtree)[,2]) -> CSTR
# # sort(inertcomp(RDAtree)[,3]) -> CSTR
# # 
# # sort(goodness(RDAtree)[,1]) -> CSTR
# # sort(goodness(RDAtree)[,2]) -> CSTR
# 
# 
# contr <- WNGmain[WNGmain$TREAT == "CONTROL",]
# CTR <- WNGmain[WNGmain$TREAT == "CONTROL",
# 
# slaSpec
# herbSpec
# ldmcSpec
#  
# # slaSub <- slaSpec[names(slaSpec) %in% names(CSTR)]
# # crtSub <- CSTR[names(CSTR) %in% names(slaSub)]
# plot(slaSub~log(crtSub[names(slaSub)]))
# slaLM <- lm(slaSub~log(crtSub[names(slaSub)]))
# summary(slaLM)
# abline(slaLM)
# 
# ldmcSub <- ldmcSpec[names(ldmcSpec) %in% names(CSTR)]
# crtSub <- CSTR[names(CSTR) %in% names(ldmcSub)]
# plot(ldmcSub~log(crtSub[names(ldmcSub)]))
# ldmcLM <- lm(ldmcSub~log(crtSub[names(ldmcSub)]))
# summary(ldmcLM)
# abline(ldmcLM)
# 
# herbSub <- herbSpec[names(herbSpec) %in% names(CSTR)]
# crtSub <- CSTR[names(CSTR) %in% names(herbSub)]
# plot(herbSub~log(crtSub[names(herbSub)]))
# herbLM <- lm(herbSub~log(crtSub[names(herbSub)]))
# summary(herbLM)
# abline(herbLM)

# Position along the partial constrained and unconstrained axes is not related with 
# trait values from control plots