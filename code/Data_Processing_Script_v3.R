# Data processing file for the Garden experiment

# Source all the datasets
source("AX_leaf_frames_analysis_data_processing.R")

# 1. Load necessary libraries
# Load all libraries
library(ggplot2)
library(vegan)
library(lme4)
library(lmerTest)
library(MASS)
library(FD)
library(xtable)
library(knitr)
library(multcomp)
library(picante)
library(devtools)
#install_github('fawda123/ggord')
#library(ggord)
#library(kableExtra)


# 2. Read all the datasets
# setwd("C:\\Users\\Piotr Szefer\\Desktop\\Work\\garden experiment\\datasets")
# treesPhylo <- read.tree("New_tree_collapsed.tre")
# #dataset <- read.csv("main_biomass_wng_ohu.csv", sep=",")
# traitData <- read.csv("tree_traits_wanang_text.csv")

# 3. Process the main dataset
# Changes in area should go here

# Replace dry mass of TREMOR with a LDMC mean:
datasetW[rownames(datasetW) == "843",]$DRY <-datasetW[rownames(datasetW) == "843",]$WET..g.*0.2896071
datasetW[rownames(datasetW) == "843",]$LDMC <- datasetW[rownames(datasetW) == "843",]$DRY/datasetW[rownames(datasetW) == "843",]$WET..g.

# Change sedge entries into shrub entries
datasetW[rownames(datasetW) == "923",]$LIFE.FORM <- "shrub"
datasetW[rownames(datasetW) == "1128",]$LIFE.FORM <- "shrub"

datasetO$SLA <- datasetO$AREA/datasetO$DRY
datasetW$SLA <- datasetW$AREA/datasetW$DRY

datasetW$HERB <- datasetW$HERB/datasetW$AREA

# logit transformation of herbivory dataset
datasetW$HERB <- log(datasetW$HERB/(1-datasetW$HERB))
datasetW$HERB[datasetW$HERB == log(0)] <- 0

dataset <- rbind(datasetO, datasetW)

# logit transformation of LDMC dataset
datasetW$LDMC <- log(datasetW$LDMC/(1-datasetW$LDMC))

# Weight shouldn't be a factor
dataset$WEIGHT <- as.numeric(as.character(dataset$WEIGHT))
datasetO$CODE <- as.character(datasetO$CODE)
datasetW$CODE <- as.character(datasetW$CODE)

datasetO$PLOT <- as.character(datasetO$PLOT)
datasetW$PLOT <- as.character(datasetW$PLOT)

datasetO$BLOCK <- as.character(datasetO$BLOCK)
datasetW$BLOCK <- as.character(datasetW$BLOCK)

# Aliases
OHUmain <- datasetO
WNGmain <- datasetW

# Changes to the SLA values
error_rows <- c("810","899","816","874","1004", 
                "1105", "1006","836","850","860",
                "1148","1160","1224","897", "947",
                "1073","936","1049","1030","1070")

WNGmain[(rownames(WNGmain) %in% error_rows), ]$SLA <- NA

WNGmain[rownames(WNGmain) == "724", ]$AREA <- 514.867
WNGmain[rownames(WNGmain) == "724", ]$SLA <- WNGmain[rownames(WNGmain) == "724", ]$AREA/WNGmain[rownames(WNGmain) == "724", ]$DRY
WNGmain[rownames(WNGmain) == "1201", ]$AREA <- 237.089
WNGmain[rownames(WNGmain) == "1201", ]$SLA <- WNGmain[rownames(WNGmain) == "1201", ]$AREA/WNGmain[rownames(WNGmain) == "1201", ]$DRY

# Replace missing values with the mean values for a given species to avoid
# over and underestimation of cwm values:

# This replaces all given species with average trait value!!!

dim(WNGmain[is.na(WNGmain$SLA), ])
WNGmain[WNGmain$SP_CODE == "CRASCR",]

# Does this even work? It should replace the NA values with mean SLA values calculted for 
# the whole experiment,
for(row in rownames(WNGmain[is.na(WNGmain$SLA), ])){
  spec <- WNGmain[rownames(WNGmain) == row,]$SP_CODE
  mean_sla <- mean(WNGmain[WNGmain$SP_CODE == spec,]$SLA, na.rm=T)
  print(spec)
  print(mean_sla)
  if(!(is.nan(mean_sla))){
    WNGmain[rownames(WNGmain) == row,]$SLA <- mean_sla
  }
}

for(row in rownames(WNGmain[is.na(WNGmain$LDMC), ])){
  spec <- WNGmain[rownames(WNGmain) == row,]$SP_CODE
  mean_ldmc <- mean(WNGmain[WNGmain$SP_CODE == spec,]$LDMC, na.rm=T)
  print(spec)
  print(mean_ldmc)
  if(!(is.nan(mean_ldmc))){
    WNGmain[rownames(WNGmain) == row,]$LDMC <- mean_ldmc
  }
}

# Check for missing values after this
# WNGmain[is.na(WNGmain$SLA), ]
# WNGmain[is.na(WNGmain$LDMC), ]

# Logged values of SLA
WNGmain$SLA <- log(WNGmain$SLA)

# Get rid of mimosa
WNGmain <- WNGmain[WNGmain$SP_CODE != "MIMODI", ]
# WNGmain[WNGmain$SP_CODE == "MIMODI", ]

# WG6P1 weevil dominateed by impereta cyllindrica!
# It also seems like the value of sla for it is too low, might be a 
# good idea to replace it with average. Entry 1156!
WNGmain[row.names(WNGmain)=="1156",]$SLA  <- 5.110306
# WG5 PREDATOR P3, high weight of herbacious plants AGERCO, SYNENO, IPOMBA
# WNGmain[WNGmain$CODE =="WG5P3",]

# Create the supplementary treatments data frame
treat <- matrix(0, nrow=length(unique(dataset$CODE)), ncol = 5)
colnames(treat) <- c("GARDEN","PLOT", "PLOT_CODE", "TREATMENT", "LOCATION")
rownames(treat) <- unique(dataset$CODE)

for (plot in unique(dataset$CODE)){
  subset <- dataset[dataset$CODE == plot,]
  treat[plot, ] <- as.matrix(subset[1, c("BLOCK","PLOT","CODE","TREAT","LOCATION")])
}

WNGtreat <- as.data.frame(treat[treat[,5] == "W", ])
# Extend with dummy treatment coding
WNGtreat$CON <- ifelse(WNGtreat$TREATMENT == "CONTROL", 1, 0)
WNGtreat$FUN <- ifelse(WNGtreat$TREATMENT == "FUNGICIDE", 1, 0)
WNGtreat$HLO <- ifelse(WNGtreat$TREATMENT == "WEEVIL25", 1, 0)
WNGtreat$HHI <- ifelse(WNGtreat$TREATMENT == "WEEVIL125", 1, 0)
WNGtreat$INS <- ifelse(WNGtreat$TREATMENT == "INSECTICIDE", 1, 0)
WNGtreat$PRE <- ifelse(WNGtreat$TREATMENT == "PREDATOR", 1, 0)
WNGtreat <- WNGtreat[order(as.character(WNGtreat$PLOT_CODE)), ]

OHUtreat <- as.data.frame(treat[treat[,5] == "O", ])
OHUtreat$CON <- ifelse(OHUtreat$TREATMENT == "CONTROL", 1, 0)
OHUtreat$FUN <- ifelse(OHUtreat$TREATMENT == "FUNGICIDE", 1, 0)
OHUtreat$HLO <- ifelse(OHUtreat$TREATMENT == "WEEVIL25", 1, 0)
OHUtreat$HHI <- ifelse(OHUtreat$TREATMENT == "WEEVIL125", 1, 0)
OHUtreat$INS <- ifelse(OHUtreat$TREATMENT == "INSECTICIDE", 1, 0)
OHUtreat$PRE <- ifelse(OHUtreat$TREATMENT == "PREDATOR", 1, 0)
OHUtreat <- OHUtreat[order(as.character(OHUtreat$PLOT_CODE)), ]


# Get the biomass data for each plot (this is, using the )
BIOSUM_W <- as.data.frame(tapply(WNGmain$WEIGHT, WNGmain$CODE, sum, na.rm = TRUE))
BIOSUM_W <- BIOSUM_W[complete.cases(BIOSUM_W),]
BIOSUM_W <- cbind(BIOSUM_W, WNGtreat)

BIOSUM_O <- as.data.frame(tapply(OHUmain$WEIGHT, OHUmain$CODE, sum, na.rm = TRUE))
BIOSUM_O <- BIOSUM_O[complete.cases(BIOSUM_O),]
BIOSUM_O <- cbind(BIOSUM_O, OHUtreat)

# This is the dataset with biomass
BIOSUM <- as.data.frame(rbind(as.matrix(BIOSUM_O),as.matrix(BIOSUM_W)))
colnames(BIOSUM) <- c("BIO", "GARDEN","PLOT", "PLOT_CODE", "TREAT", "LOCATION",
                      colnames(BIOSUM_O)[7:12])

# 2. Get the community descriptors and merge it with the BIOSUM dataset
DIV_DATA <- data.frame(GARDEN = NULL, 
                       SPEC_NO = NULL,
                       SW = NULL,
                       SIMP = NULL,
                       EVEN = NULL,
                       TREAT = NULL)

for (plot in as.character(unique(dataset$CODE))){
  subset <- dataset[dataset$CODE == plot, ]
  SPEC_NO <- length(unique(subset$SP_CODE))
  PIs <- subset$WEIGHT/sum(subset$WEIGHT)
  SW <- -sum(PIs*log(PIs))
  SIMP <- sum(PIs^2)
  EVEN <- log(SPEC_NO)
  ROW <- data.frame(GARDEN = plot,
                    SPEC_NO = SPEC_NO,
                    SW = SW,
                    SIMP = SIMP,
                    EVEN = SW/EVEN,
                    TREAT = as.character(unique(subset$TREAT)))
  DIV_DATA <- rbind(DIV_DATA, ROW)
}

DIV_DATA$BLOCK <- as.factor(substr(DIV_DATA$GARDEN, 1,3))

# HERE there is a mistake###
DESCRIPTORS <- cbind(DIV_DATA, BIOSUM)

MH_DIST <- data.frame("CONTROL"     = NULL,
                      "FUNGICIDE"   = NULL,
                      "INSECTICIDE" = NULL,
                      "PREDATOR"    = NULL,
                      "WEEVIL125"   = NULL,
                      "WEEVIL25"    = NULL)
Main <- rbind(OHUmain, WNGmain) 

for (block in unique(Main$BLOCK)){
  # Create a subset for a given block
  BCK <- Main[Main$BLOCK == block, c("TREAT", "SP_CODE", "WEIGHT")]
  
  # Get the treatments from the sample
  treatments <- sort(unique(BCK$TREAT))
  species <- sort(unique(BCK$SP_CODE))
  
  # Get the matrix of the results ready
  BCKres <- matrix(0, ncol = length(treatments), nrow = length(species))
  colnames(BCKres) <- treatments
  rownames(BCKres) <- species
  
  # loop through all the treatments and fill the columns of the results matrix
  for (treat in treatments){
    
    # Check each treatment and fill the column
    treatment_subset <- BCK[BCK$TREAT == treat,] # for loop counter for treatmnets
    
    # Loop through all the species in the given treatment and put them in the table
    for(spec in unique(treatment_subset$SP_CODE)){
      BCKres[which(rownames(BCKres) == spec), treat] <- sum(treatment_subset[treatment_subset$SP_CODE == spec, ]$WEIGHT)       # for loop counter for treatmnets
    }
  }
  # Pront the results
  MHdist <- vegdist(t(BCKres), method = "morisita", na.rm = TRUE)
  RESsub <- as.matrix(MHdist)[,1]
  MH_DIST <- rbind(MH_DIST, RESsub)
}

colnames(MH_DIST) <- treatments
rownames(MH_DIST) <- unique(Main$BLOCK)

# DESCRIPTORS dataset containing all the community descriptors
DESCRIPTORS <- DESCRIPTORS[, c(-1,-6,-7)]
DESCRIPTORS <- DESCRIPTORS[order(DESCRIPTORS$TREAT),]

# Connect MH distsance with the DESCRIPTOR DATASET
DESCRIPTORS$DIST <- as.data.frame(stack(MH_DIST))[,1]

# Fix few variables
DESCRIPTORS$BIO <- as.numeric(as.character(DESCRIPTORS$BIO))

# Biomass
DATA <- as.data.frame(DESCRIPTORS[DESCRIPTORS$TREAT == "CONTROL",c("BIO","GARDEN")])
DATA$BIO <- as.numeric(as.character(DATA$BIO))
DATA$GARDEN <- as.factor(DATA$GARDEN)
DATA$LOC <- rep(c("O","W"), c(6,6))

DATA_O <- DESCRIPTORS[DESCRIPTORS$LOCATION == "O" & DESCRIPTORS$DIST != 0, ] # 
DATA_O$TREAT <- as.factor(as.character(DATA_O$TREAT))
DATA_W <- DESCRIPTORS[DESCRIPTORS$LOCATION == "W" & DESCRIPTORS$DIST != 0, ] # 
DATA_W$TREAT <- as.factor(as.character(DATA_W$TREAT))


### Dataset for figure 1.


################
### OLD ONES ###
################

biomass <- DESCRIPTORS[DESCRIPTORS$LOCATION == "W",c("PLOT_CODE","GARDEN","PLOT","TREAT","BIO")]
biomass$TYPE <- "Biomass"
biomass$BIO <- log(biomass$BIO)
species <- DESCRIPTORS[DESCRIPTORS$LOCATION == "W",c("PLOT_CODE","GARDEN","PLOT","TREAT","SPEC_NO")]
species$TYPE <- "Species"
diversi <- DESCRIPTORS[DESCRIPTORS$LOCATION == "W",c("PLOT_CODE","GARDEN","PLOT","TREAT","SW")]
diversi$TYPE <- "Diversity"

colnames(biomass) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE","TYPE")
colnames(species) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE","TYPE")
colnames(diversi) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE","TYPE")

fig1data <- rbind(biomass, species, diversi)

## Biomass for trees, dont use cordyline terinalis and crops
#  TreeBiomass <- datasetW[datasetW$LIFE.FORM == "tree",]
TreeBiomass <- WNGmain[WNGmain$LIFE.FORM %in% c("tree", "shrub"),]

## Species dropped from the datasetW, based on the LIFE.FORMs
toDrop <- c("Manihot esculenta", "Artocarpus communis", "Carica papaya")
TreeBiomass <- TreeBiomass[!(TreeBiomass$SPEC %in% toDrop),]

## Prepare one dataset containing sum of LN(BIO), 
## STEM_NO, NO_SPEC, S-W, SIMP, sum of BASAL_A

TreeBiomass$WEIGHT <- as.numeric(TreeBiomass$WEIGHT)
TreeBiomass$BASAL_A <- as.numeric(TreeBiomass$BASAL_A)

BioTree <- tapply(TreeBiomass$WEIGHT, TreeBiomass$CODE, sum)
BasTree <- tapply(TreeBiomass$BASAL_A, TreeBiomass$CODE, sum, na.rm=T)
DivTree <- tapply(TreeBiomass$WEIGHT, TreeBiomass$CODE, function(x){-sum((x/sum(x))*log(x/sum(x)))})
SimTree <- tapply(TreeBiomass$WEIGHT, TreeBiomass$CODE, function(x){sum((x/sum(x))^2)})
SpeTree <- tapply(TreeBiomass$WEIGHT, TreeBiomass$CODE, length)
SteTree <- tapply(TreeBiomass$NO_STEMS, TreeBiomass$CODE, sum, na.rm=T)

Values <- c(log(BioTree),SpeTree,SteTree,DivTree)
Type <- rep(c("Biomass", "Richness", "Stems", "Diversity"), c(rep(36,4)))
Treatments <- rbind(WNGtreat[,c(3,4)],WNGtreat[,c(3,4)],
                    WNGtreat[,c(3,4)],WNGtreat[,c(3,4)])
TreePlotDataset <- cbind(Treatments,Values,Type)
TreePlotDataset$GARDEN <- substr(TreePlotDataset$PLOT_CODE, 1, 3)

TreeTestDataset <- cbind(WNGtreat, log(BioTree), log(BasTree+1), DivTree, 1/SimTree, 
                         SpeTree, SteTree)
colnames(TreeTestDataset) <- c("GARDEN","PLOT","PLOT_CODE","TREATMENT","LOCATION","CON","FUN",
                               "HLO","HHI","INS","PRE","logBio","logBasal","DivTree","invSimp",
                               "species", "stems")

# In the tree phylogeny I dont want to have Cordyline there (prune the tree)
TreeBiomassAsInTree <- NULL

# DATASETS for the Functional diversity
# Reset the plant names
datasetO$SP_CODE <- as.character(datasetO$SP_CODE)
WNGmain$SP_CODE <- as.character(WNGmain$SP_CODE)

# SLA is m^2 over kg (Cornelissen et al 2003) here is cm^2/g
datasetO$SLA <- datasetO$AREA/datasetO$WET..g.
traits <- tapply(datasetO$SLA, datasetO$SP_CODE, mean, na.rm=TRUE)
traits1 <- tapply(datasetO$HERB, datasetO$SP_CODE, mean, na.rm=TRUE)
traits3 <- tapply(datasetO$LDMC, datasetO$SP_CODE, mean, na.rm=TRUE)

# Repeat for wanang
traitsW <- tapply(WNGmain$SLA, WNGmain$SP_CODE, mean, na.rm = TRUE)
traitsW1 <- tapply(WNGmain$HERB, WNGmain$SP_CODE, mean, na.rm = TRUE)
traitsW3 <- tapply(WNGmain$LDMC, WNGmain$SP_CODE, mean, na.rm = TRUE)

# traitsW <- tapply(datasetW$SLA, datasetW$SP_CODE, mean, na.rm = TRUE)
# traitsW1 <- tapply(datasetW$HERB, datasetW$SP_CODE, mean, na.rm = TRUE)
# traitsW3 <- tapply(datasetW$LDMC, datasetW$SP_CODE, mean, na.rm = TRUE)


## Ordination datasets based on weight for all the species

# Source the contingency table function
source("C:\\Users\\Piotr Szefer\\Desktop\\Work\\garden experiment\\code\\contingencyTable.R")


###  FUNCTIONAL DIVERSITY 

# Trees Data
treesDBH <- WNGmain[!is.na(WNGmain$BASAL_A), c("CODE","PLOT","BLOCK", "TREAT","SPEC", "SP_CODE",
                                               "NO_STEMS", "BASAL_A", "LEAVES", "TRUNK", "WEIGHT",
                                               "HERB", "SLA", "LIFE.FORM","WATER", "LDMC")]
treesDBH <- treesDBH[treesDBH$LIFE.FORM %in% c("tree", "shrub"),]
# Removal of  Carica papaya and Manihot esculenta gives us 19 tree species
treesDBH <- treesDBH[treesDBH$SPEC != "Carica papaya", ]
treesDBH <- treesDBH[treesDBH$SPEC != "Manihot esculenta", ]
treesDBH$SPEC <- as.character(treesDBH$SPEC)
# Also change some misspeled names and assume that Premna is Premna obtusifolia
treesDBH[treesDBH$SPEC == "Ficus hispidoides", ]$SPEC <- "Ficus hispidioides"
treesDBH[treesDBH$SPEC == "Premna sp. 1", ]$SPEC <- "Premna obtusifolia"
treesDBH[treesDBH$SP_CODE == "PREMS1", ]$SP_CODE <- "PREMOB"
treesDBH[treesDBH$SPEC == "Vitex coffasus", ]$SPEC <- "Vitex cofassus"
treesDBH[treesDBH$SPEC == "Piptrus argenteus", ]$SPEC <- "Pipturus argenteus"

### COMMUNITY WEIGHTED MEANS
# For each CODE in the dataset divide weights over sum_of_weight 
WNGmain$pWeight <- 0
cwmDataTrees <- WNGmain[WNGmain$LIFE.FORM %in% c("tree", "shrub"),]

# Produces a vector of proportions (weights)
test <- tapply(WNGmain$WEIGHT, WNGmain$CODE, 
               function(X){return(X/sum(X))})
testTree <- tapply(cwmDataTrees$WEIGHT, 
                   cwmDataTrees$CODE, 
                   function(X){return(X/sum(X))})

# Stack them to obtain one vector of weights for each plot
CWMds <- stack(test) # There are 36 plots, sum = 36.
traitsDS <- WNGmain[,c("CODE","PLOT","BLOCK","TREAT","SLA","LDMC","HERB")] 

CWMdsTrees <- stack(testTree) # There are 36 plots, sum = 36.
traitsDSTrees <- cwmDataTrees[,c("CODE","PLOT","BLOCK","TREAT","SLA","LDMC","HERB")] 

# Merge two datasets based on names
dsWord <- order(traitsDS$CODE)
dsWordTrees <- order(traitsDSTrees$CODE)

CWMds <- CWMds[order(CWMds$ind),]
CWMdsTrees <- CWMdsTrees[order(CWMdsTrees$ind),]

traitsDS <- traitsDS[order(traitsDS$CODE),]
traitsDSTrees <- traitsDSTrees[order(traitsDSTrees$CODE),]

#Now they are the same no FALSE's
sum(!(traitsDS$CODE == CWMds$ind))
CWMds <- cbind(CWMds, traitsDS)

sum(!(traitsDSTrees$CODE == CWMdsTrees$ind))
CWMdsTrees <- cbind(CWMdsTrees, traitsDSTrees)

####
# Before calculating SLA and LDMC I need to know where are the missing values.

# SLA
cwm_data <- data.frame()
row <- 1
for (each_garden in unique(WNGmain$CODE)){
  subset <- WNGmain[WNGmain$CODE == each_garden,
                    c("CODE","TREAT","BLOCK","WEIGHT","SLA","LDMC","HERB")]
  #print(sum(!complete.cases(subset)))
  sub_data_sla <- subset[complete.cases(subset[,c("WEIGHT","SLA")]),c("WEIGHT","SLA")]
  #print(sub_data_sla)
  cwm_sub_data_sla <- sub_data_sla$SLA * sub_data_sla$WEIGHT/sum(sub_data_sla$WEIGHT)
  cwm_data[row, 1] <- garden <- as.character(subset[1,"BLOCK"])
  cwm_data[row, 2] <- as.character(subset[1,"TREAT"])
  cwm_data[row, 3] <- sum(cwm_sub_data_sla)
  row <- row+1
}

# plot(as.numeric(as.factor(cwm_data$V2))~cwm_data$V3, col = as.factor(cwm_data$V2),
#      pch = 19)

# LDMC
cwm_data_ldmc <- data.frame()
row <- 1
for (each_garden in unique(WNGmain$CODE)){
  subset <- WNGmain[WNGmain$CODE == each_garden,
                    c("CODE","TREAT","BLOCK","WEIGHT","SLA","LDMC","HERB")]
  #print(sum(!complete.cases(subset)))
  sub_data_ldmc <- subset[complete.cases(subset[,c("WEIGHT","LDMC")]),c("WEIGHT","LDMC")]
  #print(sub_data_sla)
  cwm_sub_data_ldmc <- sub_data_ldmc$LDMC * sub_data_ldmc$WEIGHT/sum(sub_data_ldmc$WEIGHT)
  cwm_data_ldmc[row, 1] <- as.character(subset[1,"BLOCK"])
  cwm_data_ldmc[row, 2] <- as.character(subset[1,"TREAT"])
  cwm_data_ldmc[row, 3] <- sum(cwm_sub_data_ldmc)
  row <- row+1
}

# LDMC
cwm_data_herb <- data.frame()
row <- 1
for (each_garden in unique(WNGmain$CODE)){
  subset <- WNGmain[WNGmain$CODE == each_garden,
                    c("CODE","TREAT","BLOCK","WEIGHT","SLA","LDMC","HERB")]
  sub_data_herb <- subset[complete.cases(subset[,c("WEIGHT","HERB")]),c("WEIGHT","HERB")]
  cwm_sub_data_herb <- sub_data_herb$HERB * sub_data_herb$WEIGHT/sum(sub_data_herb$WEIGHT)
  cwm_data_herb[row, 1] <- as.character(subset[1,"BLOCK"])
  cwm_data_herb[row, 2] <- as.character(subset[1,"TREAT"])
  cwm_data_herb[row, 3] <- sum(cwm_sub_data_herb)
  row <- row+1
}

cwm_data$type <- "SLA"
cwm_data_ldmc$type <- "LDMC"
cwm_data_herb$type <- "HERBIVORY"

cwm_data <- rbind(cwm_data,cwm_data_ldmc, cwm_data_herb)
names(cwm_data) <- c("garden","treat","value","type")

# plot(as.numeric(as.factor(cwm_data_ldmc$V2))~cwm_data_ldmc$V3, col = as.factor(cwm_data_ldmc$V2),
#     pch = 19)

#WG5 PREDATOR
#WNGmain[WNGmain$BLOCK =="WG5",]
# Also because high values for herbacious plants... seems like an outlier...
#WNGmain[WNGmain$CODE =="WG5P3",]
# this plot however doesn't seem unusuall...
#WNGmain[WNGmain$SP_CODE == "PIPEUM",]$LDMC
#qqnorm(WNGmain[WNGmain$SP_CODE == "MIKAMI",]$LDMC)

# Entry 1076 WG5P2 FUNGICIDE , seems like it belong to a different smaple, 
# dates dont match! see which similar plots were analysed on 7/25/16!

WNGmainTrees <- WNGmain[WNGmain$LIFE.FORM %in% c("shrub","tree"),]

#SLA TREES
cwm_data_trees <- data.frame()
row <- 1
for (each_garden in unique(WNGmainTrees$CODE)){
  subset <- WNGmainTrees[WNGmainTrees$CODE == each_garden,
                         c("CODE","TREAT","SP_CODE","BLOCK","WEIGHT","SLA","LDMC","HERB")]
  print(subset)
  sub_data_sla <- subset[complete.cases(subset[,c("WEIGHT","SLA")]),c("WEIGHT","SLA")]
  #print(sub_data_sla)
  cwm_sub_data_sla <- sub_data_sla$SLA * sub_data_sla$WEIGHT/sum(sub_data_sla$WEIGHT)
  cwm_data_trees[row, 1] <- as.character(subset[1,"BLOCK"])
  cwm_data_trees[row, 2] <- as.character(subset[1,"TREAT"])
  cwm_data_trees[row, 3] <- sum(cwm_sub_data_sla)
  row <- row+1
}

# plot(as.numeric(as.factor(cwm_data_trees$V2))~cwm_data_trees$V3, col = as.factor(cwm_data_trees$V2),
#      pch = 19)


#LDMC TREES

cwm_data_trees_ldmc <- data.frame()
row <- 1
for (each_garden in unique(WNGmainTrees$CODE)){
  subset <- WNGmainTrees[WNGmainTrees$CODE == each_garden,
                         c("CODE","TREAT","BLOCK","WEIGHT","SLA","LDMC","HERB")]
  print(subset)
  sub_data_ldmc <- subset[complete.cases(subset[,c("WEIGHT","LDMC")]),c("WEIGHT","LDMC")]
  cwm_sub_data_ldmc <- sub_data_ldmc$LDMC * sub_data_ldmc$WEIGHT/sum(sub_data_ldmc$WEIGHT)
  print(cwm_sub_data_ldmc)
  cwm_data_trees_ldmc[row, 1] <- as.character(subset[1,"BLOCK"])
  cwm_data_trees_ldmc[row, 2] <- as.character(subset[1,"TREAT"])
  cwm_data_trees_ldmc[row, 3] <- sum(cwm_sub_data_ldmc)
  row <- row+1
}


#HERB TREES

cwm_data_trees$type <- "SLA"
cwm_data_trees_ldmc$type <- "LDMC"
cwm_data_trees <- rbind(cwm_data_trees, cwm_data_trees_ldmc)
names(cwm_data_trees) <- c("garden", "treat","value","type")

####

CWMds$cwmLDMC <- CWMds$values*CWMds$LDMC
CWMds$cwmSLA <- CWMds$values*CWMds$SLA
CWMds$cwmHERB <- CWMds$values*CWMds$HERB

CWMdsTrees$cwmLDMC <- CWMdsTrees$values*CWMdsTrees$LDMC
CWMdsTrees$cwmSLA <- CWMdsTrees$values*CWMdsTrees$SLA
CWMdsTrees$cwmHERB <- CWMdsTrees$values*CWMdsTrees$HERB

# Now I need to calculate sums for the whole community
cwmLDMC <- tapply(CWMds$cwmLDMC, CWMds$CODE, sum, na.rm=TRUE)
cwmSLA <- tapply(CWMds$cwmSLA, CWMds$CODE, sum, na.rm=TRUE)
cwmHERB <- tapply(CWMds$cwmHERB, CWMds$CODE, sum, na.rm=TRUE)

# Now I need to calculate the means for trees only
cwmLDMCt <- tapply(CWMdsTrees$cwmLDMC, CWMdsTrees$CODE, sum, na.rm=TRUE)
cwmSLAt <- tapply(CWMdsTrees$cwmSLA, CWMdsTrees$CODE, sum, na.rm=TRUE)
cwmHERBt <- tapply(CWMdsTrees$cwmHERB, CWMdsTrees$CODE, sum, na.rm=TRUE)


# Get the treatment names for all
CWMtest <- cbind(WNGtreat, cwmSLA)
CWMtest$cwmLDMC <- cwmLDMC
CWMtest$cwmHERB <- cwmHERB

# Get the treatment names for trees
CWMtestTrees <- cbind(WNGtreat, cwmSLAt)
CWMtestTrees$cwmLDMCt <- cwmLDMCt
CWMtestTrees$cwmHERBt <- cwmHERBt

# Dataset for ggplot for all
SLAdat <- CWMtest[,c("TREATMENT","cwmSLA")]
colnames(SLAdat) <- c("Treat","Value")
SLAdat$Type <- "SLA"
LDMCdat <- CWMtest[,c("TREATMENT","cwmLDMC")]
colnames(LDMCdat) <- c("Treat","Value")
LDMCdat$Type <- "LDMC"
HERBdat <- CWMtest[,c("TREATMENT","cwmHERB")]
colnames(HERBdat) <- c("Treat","Value")
HERBdat$Type <- "Herbivory"
CWMplot <- rbind(SLAdat,LDMCdat,HERBdat)

# Dataset for ggplot for trees
SLAdatT <- CWMtestTrees[,c("TREATMENT","cwmSLAt")]
colnames(SLAdatT) <- c("Treat","Value")
SLAdatT$Type <- "SLA"
LDMCdatT <- CWMtestTrees[,c("TREATMENT","cwmLDMCt")]
colnames(LDMCdatT) <- c("Treat","Value")
LDMCdatT$Type <- "LDMC"
HERBdatT <- CWMtestTrees[,c("TREATMENT","cwmHERBt")]
colnames(HERBdatT) <- c("Treat","Value")
HERBdatT$Type <- "Herbivory"
CWMplotT <- rbind(SLAdatT,LDMCdatT,HERBdatT)


# I need to add traits! Check the name of the tree and add a trait for it.
traitData <- traitData[,c("SPEC", "FAMILY", "GENUS","CARB_M","NITR_M")]
traitData$FAMILY <- as.character(traitData$FAMILY)
traitData$GENUS <- as.character(traitData$GENUS)

#Add average values of SLA, HERB, WATER to the traitData dataset
# This is based on the dataset with trees above 1cm DBH for the FD's
SLA_average_trees <- tapply(treesDBH$SLA, treesDBH$SPEC, mean, na.rm=TRUE)
HERB_average_trees <- tapply(treesDBH$HERB, treesDBH$SPEC, mean, na.rm=TRUE)
WATER_average_trees <- tapply(treesDBH$WATER, treesDBH$SPEC, mean, na.rm=TRUE)
LDMC_average_trees <- tapply(treesDBH$LDMC, treesDBH$SPEC, mean, na.rm=TRUE)

# CWM for tree species in the plot. 
# I need to change names in the traitData so that they match treesDBH
# names.
traitData$SPEC <- as.character(traitData$SPEC)
rownames(traitData) <- traitData$SPEC
traitData["Melochia umbellata", ]$SPEC <- "Melochia sp. 1"
traitData["Vitex coffasus", ]$SPEC <- "Vitex cofassus"
datasetNames <- sort(unique(treesDBH$SPEC))
traitSpNames <- as.character(traitData$SPEC)

#datasetNames[!(datasetNames %in% traitSpNames)]
treesDBH$C <- traitData[treesDBH$SPEC, "CARB_M"]
treesDBH$N <- traitData[treesDBH$SPEC, "NITR_M"]



# CWM vals
ptree <- tapply(treesDBH$WEIGHT, 
                treesDBH$CODE, 
                function(X){return(X/sum(X))})
treesDBH$PROP <- stack(ptree)[,1]

treesCWM_C <- tapply(treesDBH$C*treesDBH$PROP, treesDBH$CODE, mean, na.rm=T)
treesCWM_N <- tapply(treesDBH$N*treesDBH$PROP, treesDBH$CODE, mean, na.rm=T)
treesCWM_CN <- tapply((treesDBH$C/treesDBH$N)*treesDBH$PROP, treesDBH$CODE, mean, na.rm=T)

# Faith PD
# Sample dataset: Species in columns
# Samples in rows
treesTreat <- contingencyTable2(treesDBH, "TREAT", "SPEC", "WEIGHT")
treesCode  <- contingencyTable2(treesDBH, "CODE","SPEC", "WEIGHT")

# Remove Cordyline terminalis from the dataset
treePruned <- drop.tip(treesPhylo,"'Cordylinefruticosa'")

Cnames <- colnames(treesTreat)
Cnames <- Cnames[order(order(treePruned$tip.label))]
treePruned$tip.label <- Cnames

faithPD <- pd(treesCode, treePruned)
meanPD <- mpd(treesCode, cophenetic(treePruned))
mntdPD <-mntd(treesCode, cophenetic(treePruned))

#faithPD$CODE <- rownames(faithPD)

# Put all the data together
treeDBHtest <- TreeTestDataset[,c("TREATMENT", "GARDEN")]
treeDBHtest$CODE <- rownames(treeDBHtest)

NCfPdata <- as.data.frame(cbind(treesCWM_C,treesCWM_N,treesCWM_CN,faithPD$PD,
                                meanPD, mntdPD))
NCfPdata$CODE <- rownames(NCfPdata)
colnames(NCfPdata) <- c("Carbon content CWM", "Nitrogen content CWM", 
                        "Carbon to Nitrogen content CWM","faithPD",
                        "meanPD","mntdPD")

dataToAdd <- NCfPdata[treeDBHtest$CODE, ]

## join two datasets treeDBHdata and dataToAdd
treeDBHdata <- cbind(treeDBHtest, dataToAdd)
treeDBHdata <- treeDBHdata[,-dim(treeDBHdata)[2]]
treeDBHdata[,4] <- treeDBHdata[,4]/100
treeDBHdata[,5] <- treeDBHdata[,5]/100
colnames(treeDBHdata) <- c("TREATMENT","GARDEN","CODE", 
                           "CC","NC","CtoN", "fPD", "meanPD",
                           "mntdPD")
treeDBHplot <- stack(treeDBHdata[,c("CC",
                                    "NC",
                                    "CtoN")])
treeDBHplot <- cbind(treeDBHdata[,c(1:3)],treeDBHplot)

###########################################################################################

# Functional Diversity Datasets

# Dataset 1. All species.
plantsW <- contingencyTable2(WNGmain, "CODE", "SP_CODE", "WEIGHT" )

# Dataset 2. Tree species
treesW <- contingencyTable2(WNGmain[WNGmain$LIFE.FORM %in% c("shrub","tree"),],"CODE", "SP_CODE", "WEIGHT")

# Dataset 3. Tree species > 1cm DBH, based on weight
treesWdbh <- contingencyTable2(treesDBH, "CODE", "SP_CODE", "WEIGHT")

# Traits for Wanang trees
ttreeW  <- tapply(TreeBiomass$SLA, TreeBiomass$SP_CODE, mean, na.rm = TRUE)
ttreeW1 <- tapply(TreeBiomass$HERB, TreeBiomass$SP_CODE, mean, na.rm = TRUE)
ttreeW3 <- tapply(TreeBiomass$LDMC, TreeBiomass$SP_CODE, mean, na.rm = TRUE)

# Traits for Wanang trees > 1cm DBH
ttdbhW   <- tapply(treesDBH$SLA, treesDBH$SP_CODE, mean, na.rm = TRUE)
ttdbhW1  <- tapply(treesDBH$HERB, treesDBH$SP_CODE, mean, na.rm = TRUE)
ttdbhW3  <- tapply(treesDBH$LDMC, treesDBH$SP_CODE, mean, na.rm = TRUE)

# Obtain from traitData
ttdbhWC  <- traitData[-which(rownames(traitData) == "Cordyline terminalis"),]$CARB_M
ttdbhWN  <- traitData[-which(rownames(traitData) == "Cordyline terminalis"),]$NITR_M
ttdbhWCN <- ttdbhWC/ttdbhWN

## Functional diversity for Ohu and Wanang
# a <- plantsO
aw <- plantsW
tw <- treesW
twdbh <- treesWdbh

# x <- cbind(traits, traits1, traits3)
xw <- cbind(traitsW, traitsW1, traitsW3)
xtw <- cbind(ttreeW,ttreeW1,ttreeW3) # SLA, HERB, LDMC
xtwdbh <- cbind(ttdbhW, ttdbhW1, ttdbhW3, ttdbhWC, ttdbhWN)

# Species which are absent in the communities OHU
# miss_n <- names(which(colSums(a) == 0))
# a <- as.matrix(a)
# a <- a[, -which(colnames(a) %in% miss_n)]
# x <- x[-which(rownames(x) %in% miss_n),]

# colnames(x) <- c("SLA", "HERB", "LDMC")
colnames(xw) <- c("SLA", "HERB", "LDMC")
colnames(xtw) <- c("SLA", "HERB", "LDMC")
colnames(xtwdbh) <- c("SLA", "HERB", "LDMC", "C", "N")

# Remove species with NA's in the trait dataset
# trait_names <- rownames(x[which(complete.cases(x)),])
# x <- x[complete.cases(x),]
# a <- a[, trait_names]

trait_names <- rownames(xw[which(complete.cases(xw)),])
xw <- xw[complete.cases(xw),]
aw <- aw[, trait_names]

#################################
# Introduce missing trait values:
# Here? Or in the main?
xtwdbh[!(complete.cases(xtwdbh)),"C"] <- c(mean(43.600,36.220,37.911),
                                           47.460)
xtwdbh[!(complete.cases(xtwdbh)),"N"] <- c(mean(1.840,2.160,2.122),
                                           3.180)
# Additional traits can be added here. Genus, family.
# Genus names for WNG
genuses <- substr(rownames(xw),1,4)
genusesTree <- substr(rownames(xtw),1,4)
genusesTreeDbh <- substr(rownames(xtwdbh),1,4)
# data <- as.data.frame(x)
dataW <- as.data.frame(xw)
datatW <- as.data.frame(xtw)
datatdbhW <- as.data.frame(xtwdbh)
dataW$genus <- genuses
datatW$genus <- genusesTree
datatdbhW$genus <- genusesTreeDbh

# Treatment dataset for WNG (also for the basal area dataset)

# Functional Diversity dataset
# Convert data frames into matrices
# mat_a <- matrix(0, nrow = dim(a)[1], ncol = dim(a)[2])
# colnames(mat_a) <- colnames(a)
# rownames(mat_a) <- rownames(a)
# for (row in 1:dim(a)[1]){
#   for (col in 1:dim(a)[2])
#     mat_a[row,col] <- a[row,col]
# }

dftomat <- function(aw){
  mat_aw <- matrix(0, nrow = dim(aw)[1], ncol = dim(aw)[2])
  colnames(mat_aw) <- colnames(aw)
  rownames(mat_aw) <- rownames(aw)
  for (row in 1:dim(aw)[1]){
    for (col in 1:dim(aw)[2])
      mat_aw[row,col] <- aw[row,col]
  }
  return(mat_aw)
}

mat_aw <- dftomat(aw)
mat_atw <- dftomat(tw)
mat_atwdbh <- dftomat(twdbh)

#FD <- dbFD(x[,c(1,3)],mat_a) # SLA [1] and LDMC [3].

# FDw <- dbFD(xw[,c(1,3)],mat_aw) # SLA [1] and LDMC [3].

##########
## XXX ###
##########

datatW <- datatW[,-4]
dataW <- dataW[,-4]


##########
## XXX ###
##########

FDw <- dbFD(dataW,mat_aw) # SLA [1], Herbivory, LDMC NO! genus names [3].
mat_atw_red <- mat_atw[, colnames(mat_atw)%in% rownames(datatW)]
FDtw <- dbFD(datatW,mat_atw_red) # some problem in here
FDtwdbh <- dbFD(datatdbhW,mat_atwdbh)

# Create a dataset for ggplot.

# Comparison between FD for each plot
theme <- names(FDw)[c(3,5,6,7,8)]

# Create a dataset for plotting the functional diversities plots
fdDataW <- WNGtreat[, c("PLOT_CODE","TREATMENT")] # It has 36 rows. I will make 5 graphs, so replicate 36 5 times

plotCodeW <- rep(fdDataW$PLOT_CODE, 5)
treatmentW <- rep(fdDataW$TREATMENT, 5)

dsAll <- stack(FDw[c(3,5,6,7,8)])
dsTrees <- stack(FDtw[c(3,5,6,7,8)])
dsTreesDbh <- stack(FDtwdbh[c(3,5,6,7,8)]) #WG1P6 is missing

fdDataW <- data.frame(plotCode = plotCodeW,
                      treatment = treatmentW,
                      fdVals = dsAll)
fdDataTW <- data.frame(plotCode = plotCodeW,
                       treatment = treatmentW,
                       fdVals = dsTrees)
fdDataTWdbh <- data.frame(plotCode = rownames(FDtwdbh$CWM),
                          treatment = treatmentW[-(which(names(treatmentW)=="WG1P6"))],
                          fdVals = dsTreesDbh)

########################################
### ORDERING OF fdDataW and fdDataTW ###
########################################

# Change names of treatments
fdDataW$treatment <- as.character(fdDataW$treatment)
fdDataW[fdDataW$treatment == "PREDATOR",]$treatment <- "P"
fdDataW[fdDataW$treatment == "INSECTICIDE",]$treatment <- "I"
fdDataW[fdDataW$treatment == "FUNGICIDE",]$treatment <- "F"
fdDataW[fdDataW$treatment == "WEEVIL25",]$treatment <- "H1"
fdDataW[fdDataW$treatment == "WEEVIL125",]$treatment <- "H2"
fdDataW[fdDataW$treatment == "CONTROL",]$treatment <- "C"

treat_ord <- c("F","I","C","P","H1","H2")
#treat_ord <- c("C","F","I","P","H1","H2")
type_ord <- c("FRic","FEve","FDiv","FDis","RaoQ")

fdDataW <- fdDataW[order(match(as.character(fdDataW$treatment),treat_ord),
                             match(as.character(fdDataW$fdVals.ind),type_ord)),]

# It hacs to be sorted like in the graph, so that the colours work
fdDataW$treatment <- factor(fdDataW$treatment, levels = treat_ord)
fdDataW$fdVals.ind <- factor(fdDataW$fdVals.ind, levels = type_ord)
fdDataW <- fdDataW[order(match(as.character(fdDataW$fdVals.ind),
                               type_ord)), ]

###### fdDataTW

# Change names of treatments
fdDataTW$treatment <- as.character(fdDataTW$treatment)
fdDataTW[fdDataTW$treatment == "PREDATOR",]$treatment <- "P"
fdDataTW[fdDataTW$treatment == "INSECTICIDE",]$treatment <- "I"
fdDataTW[fdDataTW$treatment == "FUNGICIDE",]$treatment <- "F"
fdDataTW[fdDataTW$treatment == "WEEVIL25",]$treatment <- "H1"
fdDataTW[fdDataTW$treatment == "WEEVIL125",]$treatment <- "H2"
fdDataTW[fdDataTW$treatment == "CONTROL",]$treatment <- "C"

treat_ord <- c("F","I","C","P","H1","H2")
#treat_ord <- c("C","F","I","P","H1","H2")
type_ord <- c("FRic","FEve","FDiv","FDis","RaoQ")

fdDataTW <- fdDataTW[order(match(as.character(fdDataTW$treatment),treat_ord),
                         match(as.character(fdDataTW$fdVals.ind),type_ord)),]

# It hacs to be sorted like in the graph, so that the colours work
fdDataTW$treatment <- factor(fdDataTW$treatment, levels = treat_ord)
fdDataTW$fdVals.ind <- factor(fdDataTW$fdVals.ind, levels = type_ord)
fdDataTW <- fdDataTW[order(match(as.character(fdDataTW$fdVals.ind),
                               type_ord)), ]

######################################################################################


# Data constructed to perform tests for Wanang data set
fdDataTest <- data.frame(Plot = WNGtreat[, c("PLOT_CODE")], Block = substr(WNGtreat[, c("PLOT_CODE")], 1,3), Treat = WNGtreat[, c("TREATMENT")],
                         FRic = FDw$FRic, FEve = FDw$FEve, FDiv = FDw$FDiv, FDis = FDw$FDis,
                         RaoQ = FDw$RaoQ)
fdDataTwTest <- data.frame(Plot = WNGtreat[, c("PLOT_CODE")], Block = substr(WNGtreat[, c("PLOT_CODE")], 1,3), Treat = WNGtreat[, c("TREATMENT")],
                           FRic = FDtw$FRic, FEve = FDtw$FEve, FDiv = FDtw$FDiv, FDis = FDtw$FDis,
                           RaoQ = FDtw$RaoQ)
fdDataTwDbhTest <- data.frame(Plot = rownames(FDtwdbh$CWM), 
                              Block = substr(rownames(FDtwdbh$CWM), 1,3), 
                              Treat = WNGtreat[-(which(rownames(WNGtreat)=="WG1P6")),"TREATMENT"],
                              FRic = FDtwdbh$FRic, FEve = FDtwdbh$FEve, 
                              FDiv = FDtwdbh$FDiv, FDis = FDtwdbh$FDis,
                              RaoQ = FDtwdbh$RaoQ)

### Dataset for figure 1.

################
### NEW ONES ###
################

AllTestData <- DESCRIPTORS[substr(DESCRIPTORS$GARDEN,1,2) == "WG", ]
AllTestData$EVEN <- AllTestData$SW/log(AllTestData$SPEC_NO)

# Tree community TreePlotDataset
TreeTestDataset$Evenness <- TreeTestDataset$DivTree/log(TreeTestDataset$species)

# Now make a dataset that can be used for the panels
d1  <- AllTestData[,c("PLOT_CODE","GARDEN","PLOT","TREAT","BIO")]
d2  <- AllTestData[,c("PLOT_CODE","GARDEN","PLOT","TREAT","SW")]
d3  <- AllTestData[,c("PLOT_CODE","GARDEN","PLOT","TREAT","SPEC_NO")]
d4  <- AllTestData[,c("PLOT_CODE","GARDEN","PLOT","TREAT","SIMP")]
d5  <- AllTestData[,c("PLOT_CODE","GARDEN","PLOT","TREAT","EVEN")]
dt1 <- TreeTestDataset[,c("PLOT_CODE","GARDEN","PLOT","TREATMENT","logBio")]
dt2 <- TreeTestDataset[,c("PLOT_CODE","GARDEN","PLOT","TREATMENT","DivTree")]
dt3 <- TreeTestDataset[,c("PLOT_CODE","GARDEN","PLOT","TREATMENT","species")]
dt4 <- TreeTestDataset[,c("PLOT_CODE","GARDEN","PLOT","TREATMENT","invSimp")]
dt5 <- TreeTestDataset[,c("PLOT_CODE","GARDEN","PLOT","TREATMENT","Evenness")]
dt6 <- TreeTestDataset[,c("PLOT_CODE","GARDEN","PLOT","TREATMENT","stems")]

colnames(d1) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
d1$VALUE <- log(d1$VALUE)
colnames(d2) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
colnames(d3) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
colnames(d4) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
colnames(d5) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
colnames(dt1) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
colnames(dt2) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
colnames(dt3) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
colnames(dt4) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
dt4$VALUE <- 1/dt4$VALUE
colnames(dt5) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")
colnames(dt6) <- c("PLOT_CODE","GARDEN","PLOT","TREAT","VALUE")


panelData <- rbind(d1,d2,d3,d4,d5,dt1,dt2,dt3,dt4,dt5,dt6)
panelData$TYPE <- rep(c("LN(Biomass)", "Shannon's H", "Richness", "Simpson's D",
                        "Evenness",
                        "LN(Biomass) (Woody plants)", "Shannon's H (Woody plants)", 
                        "Richness (Woody plants)", "Simpson's D (Woody plants)",
                        "Evenness (Woody plants)", "No. of stems (Woody plants)"), 
                      each=36)

# Change names of treatments
panelData$TREAT <- as.character(panelData$TREAT)
panelData[panelData$TREAT == "PREDATOR",]$TREAT <- "P"
panelData[panelData$TREAT == "INSECTICIDE",]$TREAT <- "I"
panelData[panelData$TREAT == "FUNGICIDE",]$TREAT <- "F"
panelData[panelData$TREAT == "WEEVIL25",]$TREAT <- "H1"
panelData[panelData$TREAT == "WEEVIL125",]$TREAT <- "H2"
panelData[panelData$TREAT == "CONTROL",]$TREAT <- "C"

################
### ORDERING ###
################

treat_ord <- c("F","I","C","P","H1","H2")
#treat_ord <- c("C","F","I","P","H1","H2")
type_ord <- c("LN(Biomass)","Shannon's H","Richness","Simpson's D","Evenness",
              "LN(Biomass) (Woody plants)", "Shannon's H (Woody plants)",
              "Richness (Woody plants)", "Simpson's D (Woody plants)", "Evenness (Woody plants)",
              "No. of stems (Woody plants)")

panelData <- panelData[order(match(as.character(panelData$TREAT),treat_ord),
                           match(as.character(panelData$TYPE),type_ord)),]

# It hacs to be sorted like in the graph, so that the colours work
panelData$TREAT <- factor(panelData$TREAT, levels = treat_ord)
panelData$TYPE <- factor(panelData$TYPE, levels = type_ord)
panelData <- panelData[order(match(as.character(panelData$TYPE),type_ord)), ]


# Community weighted mean data for the panel
CWMplot
CWMplotT$Type <- rep(c("SLA (Woody plants)", "LDMC (Woody plants)",
                       "Herbivory (Woody plants)"), each=36)
panelCWM <- rbind(CWMplot, CWMplotT)

panelCWM$GARDEN <- substr(rownames(panelCWM),1,3)
panelCWM$cols <- "gray40"
# Change names of treatments
panelCWM$Treat <- as.character(panelCWM$Treat)
panelCWM[panelCWM$Treat == "PREDATOR",]$Treat <- "P"
panelCWM[panelCWM$Treat == "INSECTICIDE",]$Treat <- "I"
panelCWM[panelCWM$Treat == "FUNGICIDE",]$Treat <- "F"
panelCWM[panelCWM$Treat == "WEEVIL25",]$Treat <- "H1"
panelCWM[panelCWM$Treat == "WEEVIL125",]$Treat <- "H2"
panelCWM[panelCWM$Treat == "CONTROL",]$Treat <- "C"

####################
### SORTING DTFR ###
####################

treat_ord <- c("F","I","C","P","H1","H2")
#treat_ord <- c("C","F","I","P","H1","H2")
type_ord <- c("SLA", "LDMC", "Herbivory",
              "SLA (Woody plants)", "LDMC (Woody plants)",
              "Herbivory (Woody plants)")
panelCWM <- panelCWM[order(match(as.character(panelCWM$Treat),treat_ord),
                           match(as.character(panelCWM$Type),type_ord)),]

# Change the order of the factors
panelCWM$Treat <- factor(panelCWM$Treat, levels = unique(panelCWM$Treat))
panelCWM$Type <- factor(panelCWM$Type,
                        levels = unique(panelCWM$Type))
# Sort like in th eplot
panelCWM <- panelCWM[order(match(as.character(panelCWM$Type),type_ord)),]

setwd("C:/Users/Piotr Szefer/Desktop/Work/garden experiment/datasets/randomizations/trees_assembly")
dataRC_facets_plot <- as.data.frame(read.table("randomizationsTrees.txt"))
dataRC_facets_plot$treat <- as.character(dataRC_facets_plot$treat)
dataRC_facets_plot[dataRC_facets_plot$treat == "PREDATOR",]$treat <- "P"
dataRC_facets_plot[dataRC_facets_plot$treat == "INSECTICIDE",]$treat <- "I"
dataRC_facets_plot[dataRC_facets_plot$treat == "FUNGICIDE",]$treat <- "F"
dataRC_facets_plot[dataRC_facets_plot$treat == "WEEVIL25",]$treat <- "H1"
dataRC_facets_plot[dataRC_facets_plot$treat == "WEEVIL125",]$treat <- "H2"
dataRC_facets_plot[dataRC_facets_plot$treat == "CONTROL",]$treat <- "C"

############################
## Only herbaceous plants ##
############################

