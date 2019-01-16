# Randomizations

# For each randomization we constrained our dataset only to control plots and one 
# treatment plot, so that experimental communities will not be overrepresented in 
# the randomizations. Communitties for a given experimental plot were randomly 
# assembled in a following way: (1) For a given plot we restricted number of 
# species and total biomass in randomly assembled community to its empirical ones. 
# (2) Probablity of each species to be sampled to a randomized community was 
# calculated based on occurence probability in all plots for a given pair of 
# treatments. obtained from all 12 plots for a given pair of control and treatment 
# blocks.  In comparison with Stegen et al. (2013) we used individuals described by their
# biomass. (3) Cumulative biomass for a given plot was an empirical bimass 
# for that plot. Units of biomass were randomly assigned into species depending 
# on its relative biomass. (4) Probability of increasing biomas of a given speices 
# wa proportional to the total biomass of that speciess accros all 12 plots. 


# 1. Setup 

####
# Set the directory containing the data here  <<<<<<<<<<
####

# Datasets: plant_mat.txt, plant_mat_all.txt, WNGtreat.txt

path <- "C:\\Users\\Piotr Szefer\\Desktop\\Work\\garden experiment\\datasets\\randomizations"
setwd(path)
rands <- 99   # sets the number of iterations for the randomizations
accu <- 0.001 # sets the accuracy at which biomass increases during assembly
## 1a. Read the data
### Community matrix
plant_mat <- read.table("plant_mat.txt")
plant_mat_all <- read.table("plant_mat_all.txt")

### All possible paorwise comparisons of plots
#all_pairs <- read.table("all_pairs_trees.txt")
#colnames(all_pairs) <- c("P1","P2")
### Treatments
WNGtreat <- read.table("WNGtreat.txt")

####
# Install necessary packages  <<<<<<<<<<
####

#install.packages(c("vegan","stringr","dplyr", "ggplot2"))

## 1b. Load packages

library(vegan)
library(stringr)
library(dplyr)
library(ggplot2)

################################################################################
############################   FUNCTIONS #######################################
################################################################################

### Assemble function
assemble <- function(spec_patterns, bio_patterns, plot_name, plant_names,
                     rel_probs, rel_plant_biomass, accuracy = 0.001){
  
  # Spec patterns are the number of species for each of the plot
  # plot1 - character string with the name of first 
  # bio patterns biomass rounded to 3 decimal points
  # plant_names names of the plant in the whole community
  # rel_probs is numer of incidences divided by all the plots
  # rel_plant biomass ir relative biomass of each speceis within the community
  
  spec <- spec_patterns[names(spec_patterns) == plot_name]
  bio <- bio_patterns[names(bio_patterns) == plot_name]
  
  # Sample given number of species
  S1 <- sample(plant_names,spec,replace=FALSE,prob = rel_probs)
  
  # Probabilities of observing one unit of biomass for a given community
  sub_rel_plant_biomass <- rel_plant_biomass[S1]/sum(rel_plant_biomass[S1])
  
  # Assembling the plot
  initS1 <- S1 # initiate minimal biomass
  iterS1 <- sample(S1, ((bio/accuracy)-spec), TRUE, 
                   prob=sub_rel_plant_biomass)
  sampleS1 <- c(initS1, iterS1)
  randCom1 <- table(sampleS1)*accuracy
  return(list(COM = randCom1))
}

### Contingency table function
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

### Raup-Crick function

RCrands <- function(p1, p2, plant_mat,
                    empiricalBCval,
                    rands=10, accu=0.1){
  
  # 1 Calculate parameters for randomizations
  spec_patterns <- rowSums(plant_mat > 0)
  bio_patterns <- round(rowSums(plant_mat),-log(accu)/log(10))
  rel_biomass <- colSums(plant_mat)/sum(colSums(plant_mat))
  probs <- round(colSums(plant_mat>0)/dim(plant_mat)[1],3)
  rel_probs <- probs/sum(probs) # relative probabilities
  
  # 2 Randomize
  
  iterations <- rands
  bcVals <- c()
  
  for(iter in 1:iterations){
    P1 <- assemble(spec_patterns, bio_patterns, p1, names(rel_probs),
                   rel_probs, rel_biomass, accuracy = accu)
    P2 <- assemble(spec_patterns, bio_patterns, p2, names(rel_probs),
                   rel_probs, rel_biomass, accuracy = accu)
    
    # Put both into one dataset to calculate randomizedBC
    p1df <- as.data.frame(P1$COM)
    p1df$com <- "p1"
    p2df <- as.data.frame(P2$COM)
    p2df$com <- "p2"
    
    df <- rbind(p1df, p2df)
    colnames(df) <- c("spec", "bio", "com")
    p1p2 <- contingencyTable2(df,"com","spec","bio")
    
    # BC values for all randomizations
    bcVals <- c(bcVals, vegdist(p1p2, method = "bray"))
  }
  
  # Proportion of BC distance smaller or equal to the empirical one 
  # (add to pairs)
  
  randProp <- sum(bcVals<=empiricalBCval)/iterations
  
  ## 2b. Calculate the R-C index
  rcVal <- (randProp - 0.5)*2
  
  return(list(randBC=bcVals, RC=rcVal))
}

# Helper function for obtaining empirical values from BC matrix
getEmpiricalBC <- function(p1,p2,empiricalBC){
  # p1 and p2 as characters and empirical BC, vedgdist
  matBC <- as.matrix(empiricalBC)
  val <- matBC[rownames(matBC) == p1,
               colnames(matBC) == p2]
  return(val)
}


#################################################################################
######### RANDOMIZATIONS ########################################################
#################################################################################

# 2. Pairwise comparisons

# Create directory
dir.create(file.path(path, "trees_assembly"))
setwd("./trees_assembly")

## Onle these interests me. Subset the plant_mat matrix so that it contains only
## Controll and the other treatment. I need to split the dataset into smaller
## datasets

trtPairs <- expand.grid(unique(WNGtreat$TREATMENT),unique(WNGtreat$TREATMENT))
colnames(trtPairs) <- c("trt1", "trt2")
trtPairs <- trtPairs[trtPairs$trt1 != trtPairs$trt2, ]
trtPairs$trt1 <- as.character(trtPairs$trt1)
trtPairs$trt2 <- as.character(trtPairs$trt2)
trtPairs <- trtPairs[!duplicated(t(apply(trtPairs, 1, sort))),]
WNGtreat$TREATMENT <- as.character(WNGtreat$TREATMENT)
WNGtreat$PLOT_CODE <- as.character(WNGtreat$PLOT_CODE)
# Reduce treatment pairs to onle these which have CONTROL in it
trtPairs <- trtPairs[trtPairs$trt1 == "CONTROL" | trtPairs$trt2 == "CONTROL", ]
lowerPart <- cbind(trtPairs[3:5,2],
                   trtPairs[3:5,1])
colnames(lowerPart) <- colnames(trtPairs)
trtPairs <- rbind(trtPairs[1:2,1:2], 
                  lowerPart)


############################
## SUBSETTING STARTS HERE ##
############################

# comb = 1 # for testing procedures inside the function

randomization <- function(plant_mat){
  for (comb in 1:dim(trtPairs)[1]){
    # use row from trtPairs
    rn <- rownames(WNGtreat[(WNGtreat$TREATMENT %in% c(trtPairs$trt1[comb],
                                                       trtPairs$trt2[comb])),])
    print(c(trtPairs$trt1[comb],trtPairs$trt2[comb]))
    
    # Subset the community matrix
    sub_plant_mat <- plant_mat[rownames(plant_mat) %in% rn, ]
    
    # Should I give small probability to the rest of the plants to get to
    # the plot? To ones which are either not in control or not in the treatment?
    # Does this seem like a missing information?
    sub_plant_mat <- sub_plant_mat[, colSums(sub_plant_mat) != 0]
    
    # Subset the trestment data
    sub_treat <- WNGtreat[rownames(WNGtreat) %in% rn,]
    sub_treat$TREATMENT <- as.character(sub_treat$TREATMENT)
    sub_treat$PLOT_CODE <- as.character(sub_treat$PLOT_CODE)
    sub_treat$GARDEN <- as.character(sub_treat$GARDEN)
    
    # Make a matrix with colnames and rownames as plots
    RCmat <- matrix(ncol = dim(sub_plant_mat)[1], nrow=dim(sub_plant_mat))
    colnames(RCmat) <- rownames(sub_plant_mat)
    rownames(RCmat) <- rownames(sub_plant_mat)
    
    # Get all the pairs which suppsed to be filled in to RCmat, lower triangle
    indices <- which(lower.tri(RCmat), arr.ind = T)
    sub_all_pairs <- data.frame(P1 = rownames(RCmat)[indices[,1]],
                                P2 = colnames(RCmat)[indices[,2]])
    
    sub_all_pairs$P1 <- as.character(rownames(RCmat)[indices[,1]])
    sub_all_pairs$P2 <- as.character(colnames(RCmat)[indices[,2]])
    
    # BRAY MATRIX for the subset of two treatments
    empiricalBC <- vegdist(sub_plant_mat)
    
    # Now fill this matrix with RC indices
    # Loop through all the plots in sub_all_pairs and calculate RCvals for them
    
    # comp = 2
    
    for(comp in 1:dim(sub_all_pairs)[1]){
      
      p1 <- sub_all_pairs[comp,]$P1
      p2 <- sub_all_pairs[comp,]$P2
      
      # See which pair is currently randomized
      print(c(p1,p2))
      
      # Get empirical values of BC for the subsetted matrix
      empBCval <- getEmpiricalBC(sub_all_pairs$P1[comp], 
                                 sub_all_pairs$P2[comp],
                                 empiricalBC)
      
      #######################################
      # Raup Crick calculations with result #
      #######################################
      
      # Assembly communitites for a given pair
      RCresult <- RCrands(p1,p2,sub_plant_mat,
                          empBCval, rands = rands, 
                          accu = accu)
      
      #######################################
      # Put that result in the matrix RCmat #
      #######################################
      
      RCmat[rownames(RCmat) == p1, colnames(RCmat) == p2] <- RCresult$RC
    }
    # Save the comparison resulting matrix
    
    write.table(RCmat, paste("mat_",trtPairs$trt1[comb], "_",
                             trtPairs$trt2[comb],
                             ".txt",sep=""))
  }
}

randomization(plant_mat)

## Comparisons for the whole community

setwd("..")
dir.create(file.path(path, "all_assembly"))
setwd("./all_assembly")

randomization(plant_mat_all)

##################################################################################
##############  PROCESS RESULTS OF THE RANDOMIZATION #############################
##################################################################################

plotMeMyGraphs <- function(path){
  
  setwd(path)
  
  # Check all the txt files with a word CONTROL in them (RC matrices)
  file_names <- list.files()[str_detect(list.files(), "CONTROL")]
  
  # Data for facets and to store results, mean, lose, upse, type
  dataRC_facets <- data.frame()
  
  # For each file name calculate mean values of inter treatmet comparisons
  for(file in 1:length(file_names)) {
    rcmatrix <- read.table(file_names[file])
    
    # Fill the whole matrix with values from th elower triangle
    rcmatrix[upper.tri(rcmatrix)] <- rcmatrix[lower.tri(rcmatrix)]
    
    # Subset the treatment dataset only to these plots
    sub_treat <- WNGtreat[WNGtreat$PLOT_CODE %in% rownames(rcmatrix), ]
    
    # Now for a given plot of a specific treatment, take mean RCval of this
    # plot with all the others from the same treatment
    resRC <- data.frame()
    
    for (row in 1:dim(sub_treat)[1]){
      plt <- sub_treat$PLOT_CODE[row]
      trt <- sub_treat$TREATMENT[row]
      grd <- sub_treat$GARDEN[row]
      
      # Print current row
      print(c(plt,trt,grd))
      
      # Other plots belonging to the same treatment
      trtPlots <- as.character(sub_treat[sub_treat$TREATMENT == trt, ]$PLOT_CODE)
      
      print(trtPlots)
      
      mRC <- mean(as.numeric(rcmatrix[rownames(rcmatrix) == plt,
                                      (colnames(rcmatrix) %in% trtPlots)]), na.rm=T)
      
      print(mRC)
      
      # Add reults to the resRC data
      rowRC <- data.frame(plot=plt,
                          treat = trt,
                          gard = grd,
                          meanRC = mRC)
      resRC <- rbind(resRC, rowRC)
    }
    resRC$type <- paste(file_names[file])
    dataRC_facets <- rbind(dataRC_facets,resRC)
  }
  
  
  ##########################
  # 3b. ggplot facet plots #
  ##########################
  
  dataRC_facets
  dataRC_facets_plot <- data.frame()
  
  for (comparison in unique(dataRC_facets$type)){
    subDat <- dataRC_facets[dataRC_facets$type == comparison,]
    subDat$treat <- as.character(subDat$treat)
    subDat$gard <- as.character(subDat$gard)
    subDat$plot <- as.character(subDat$plot)
    
    meanRc <- tapply(subDat$meanRC, subDat$treat, mean)
    Vars   <- tapply(subDat$meanRC, subDat$treat, var)
    upSE   <- meanRc + sqrt(Vars)
    loSE   <- meanRc - sqrt(Vars)
    
    subRes <- data.frame(meanRC = meanRc, lose = loSE, upse = upSE,
                         treat = names(meanRc),
                         type = c(comparison, comparison))
    
    dataRC_facets_plot <- rbind(dataRC_facets_plot,
                                subRes)
    
  }
  
  dataRC_facets_plot$treat <- as.character(dataRC_facets_plot$treat)
  dataRC_facets_plot$type <- as.character(dataRC_facets_plot$type)
  dataRC_facets_plot$type <- gsub("mat_", "", dataRC_facets_plot$type)
  dataRC_facets_plot$type <- gsub(".txt", "", dataRC_facets_plot$type)
  dataRC_facets_plot$type <- gsub("_", " vs ", dataRC_facets_plot$type)
  
  pl <- ggplot() +
    geom_errorbar(data=dataRC_facets_plot,
                  mapping=aes(x=treat, ymin=lose, ymax=upse), width=0.2,
                  size=1, color="gray60") + ylab("Within treatment dissimilarity (RC)") +
    geom_hline(yintercept = c(0), lty=4, lwd = 1.1, col="grey60") +
    geom_point(data=dataRC_facets_plot, mapping=aes(x=treat, y=meanRC), size=4, shape=21,
               fill="black")+
    ylim(-1,1) + geom_hline(yintercept = c(1, -1), lty=2, col="grey60") +
    facet_wrap(~type, scales="free", drop=TRUE) + theme_bw()+
    theme(axis.text.x=element_text(angle=90, hjust=1))
  
  #windows(800,600)
  
  pl
  
  # Test for significance
  
  dataRC_facets$ftreat <- as.factor(dataRC_facets$treat)
  for (comparison in unique(dataRC_facets$type)){
    subDat <- dataRC_facets[dataRC_facets$type == comparison,]
    l0 <- lmer(meanRC~1+(1|gard),subDat)
    l1<- lmer(meanRC~ftreat+(1|gard),subDat)
    print(comparison)
    print(anova(l0,l1))
  }
  return(list(plot = pl,df = as.data.frame(dataRC_facets_plot)))
}


###############################
# Read the data from the file #
###############################

setwd("..")
setwd("all_assembly")
path1 <- getwd()
setwd("..")
setwd("trees_assembly")
path2 <- getwd() 

plot1 <- plotMeMyGraphs(path1)
windows(800,600)
plot1$plot

plot2 <- plotMeMyGraphs(path2)
windows(800,600)
plot2$plot

#write.table(plot2$df, "randomizationsTrees.txt")


pdf()
plot1
plot2
dev.off()