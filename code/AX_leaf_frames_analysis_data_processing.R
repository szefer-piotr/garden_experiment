# Ohu leaf frames processing file. It takes the raw data with all the leaves,
# obtained with the use of SUPL_gettog_final.R script. The task is to match teh names
# and create a new columns for the main dataset


## Read the dataset and drop the basal area columns
### Read the main dataset
dataset <- read.csv("datasets/main_biomass_wng_ohu.csv", sep=",", dec=".", header=T)
dataset <- dataset[, 1:23]

# Add location variable and fix the classes
dataset$LOCATION <- substr(dataset$CODE, 1,1)

## Weight shouldn't be a factor
dataset$WEIGHT <- as.numeric(as.character(dataset$WEIGHT))
dataset$DRY <- as.numeric(as.character(dataset$DRY))

## Read the leaf frames dataset
leafOhu <- read.csv("datasets/ohu_leaf_frames.csv", sep="\ ")
leafWng <- read.csv("datasets/wng_leaf_frames.csv", sep="\ ") # Created with getOg script

## Unify the names
leafWng$Label <- as.character(leafWng$Label)
leafWng$Label[leafWng$Label == "W11G2P6MELAMU_P.JPG"] <- "W1G2P6MELAMU_P.JPG"
leafWng$Label[leafWng$Label == "W1C4D6MACATA_P.JPG"] <-  "W1G4P6MACATA_P.JPG"
leafWng$Label[leafWng$Label == "W1G1D2PASPALUMSP1_P.JPG"] <- "W1G1P2PASPALUMSP1_P.JPG"
leafWng$Label[leafWng$Label == "W1C4P6CUCUBITASP1_P.JPG"] <- "W1G4P6CUCUBITASP1_P.JPG"
### Replace "W1" with W
leafWng$Label <- gsub("W1", "W", leafWng$Label)
leafWng$Label <- gsub("WI", "W", leafWng$Label)

#############################################################################
#############################################################################

# Using only leaf for wanang
# I can repeat the whole procedure for Ohu as well
leaf <- leafWng

#############################################################################
#############################################################################

# Area y is the real area.
# Herbivory is herbivory. I can get mean percentage tissue loss or ratio.
# Get the mean values for herbivory and sum of the weight and area (Area.y).

#meanHerb <- tapply(leafOhu$Herbivory, leafOhu$Label, mean, na.rm = TRUE) 
meanHerb <- tapply(leaf$Herbivory, leaf$Label, mean, na.rm = TRUE) 

# Hypothesized area without herbivory damage
#meanAreaFull <- tapply(leafOhu$Area.x, leafOhu$Label, sum, na.rm = TRUE) 
meanAreaFull <- tapply(leaf$Area.x, leaf$Label, sum, na.rm = TRUE) 

# Real area from the leaf frames
#meanAreaReal <- tapply(leafOhu$Area.y, leafOhu$Label, sum, na.rm = TRUE) 
meanAreaReal <- tapply(leaf$Area.y, leaf$Label, sum, na.rm = TRUE) 

# Joined big dataset for Ohu and Wanang
means <- as.data.frame(cbind(meanHerb, meanAreaFull, meanAreaReal))
means$filename <- rownames(means)
means$filename <- substr(means$filename, 1, nchar(means$filename)-6) # Drop the _P.jpg part
means$plantcode <- substr(means$filename, 6, nchar(means$filename)) # Extract plant names from photos
means$garden <- substr(means$filename,1, 3) 


# Two first rows have 0 instead of O.
#means$garden[c(1,2)] <- "OG2"
means$plot <- substr(means$filename, 4, 5)
means$code <- paste(means$garden, means$plot, sep="")

# First clean the data and then merge them
datasetW <- dataset[dataset$LOCATION == "W",]

# From the original dataset, what are the species names
specNamesWanang <- as.matrix(sort(unique(datasetW$SPEC)))
specCodesWanang <- as.matrix(sort(unique(datasetW$SP_CODE)))

spComp <- cbind(specCodesWanang,rbind(specNamesWanang,1))
datasetW$SP_CODE <- as.character(datasetW$SP_CODE)

# This is Dracantomeolon lanceolatum, DRACLA
datasetW[datasetW$SP_CODE == "DRACDA",]$SP_CODE <- "DRACLA"

# This is Dioscorea sp. 1, should be DIOSS1
datasetW[datasetW$SP_CODE == "DIOSAL",]$SP_CODE <- "DIOSS1"

datasetW[datasetW$CODE == "WG4P6" & datasetW$SP_CODE == "MACAAL",]$BLOCK <- "WG5"
datasetW[datasetW$CODE == "WG4P6" & datasetW$SP_CODE == "MACAAL",]$PLOT <- "P2"
datasetW[datasetW$CODE == "WG4P6" & datasetW$SP_CODE == "MACAAL",]$CODE <- "WG5P2"


# Compare names and cosed again
specNamesWanang <- as.matrix(sort(unique(datasetW$SPEC)))
specCodesWanang <- as.matrix(sort(unique(datasetW$SP_CODE)))
spComp <- cbind(specCodesWanang,rbind(specNamesWanang))

# Calculate the LDMC
datasetW$LDMC <- datasetW$DRY/datasetW$WET..g.

# sum(is.na(datasetW$LDMC))/length(datasetW$LDMC)  # 3.5 % missing data.

# Get the names from the leaf frame data
#as.matrix(sort(unique(means$plantcode)))
means$plantcode <- as.character(means$plantcode)

##################################
## From the manual changes file ##
##################################

means[means$plantcode == "HOMA",]$plantcode <- "HOMANO"
means[means$plantcode == "5PIPTAR",]$plot <- "P5" 
means[rownames(means) == "WG1P6FUPHHI_P.JPG",]$plantcode <- "EUPHHI"
means[rownames(means) == "WG2P2MELAMU_P.JPG",]$plot <- "P1"
means[rownames(means) == "WG2P2MELAMU_P.JPG",]$code <- "WG2P1"
means[rownames(means) == "WG2P5SIDARH_P.JPG",]$plantcode <- "MELAMU"
# means[rownames(means) == "WG4P1FICUCO_P.JPG",]$plantcode <- "FICUCO"
# means[rownames(means) == "WG4P2FICUCO_P.JPG",]$plantcode <- "FICUCO"

means[rownames(means) == "WGP2AGERCO_P.JPG", ]$code   <- "WG6P2"
means[rownames(means) == "WGP2AGERCO_P.JPG", ]$plot   <- "P2"
means[rownames(means) == "WGP2AGERCO_P.JPG", ]$garden <- "WG6"

means[rownames(means) == "WGP3CALOMU_P.JPG", ]$code   <- "WG1P3"
means[rownames(means) == "WGP3CALOMU_P.JPG", ]$plot   <- "P3"
means[rownames(means) == "WGP3CALOMU_P.JPG", ]$garden <- "WG1"

# Now change the mC file as well!!! # it needs to leave all of the above
# unchanged!
# But what should I change???
# 1. I have removed the row where FUPHHI was changet into FICUHI
# 2. after i have changed the input files in my fie there should be SNENO removed
#    form the mC dataset
# 3. Ficus congesta keeps being changed to FICUCP


# Change the plantcode in means dataset according to "meansCorrect.txt"
mC <- read.table("datasets/meansCorrect.csv", header = FALSE, sep="\t")
mC[,1] <- as.character(mC[,1])
mC[,2] <- as.character(mC[,2])

#mC[136,2] <- "TREMS1"
##################################
##################################
##################################

### For loop to correct the codes in 'means' dataset
means$done <- FALSE

#means[(means$plantcode == "AGERCO") & !means$done,]

for (code in mC[,1]){
  #print(code)
  index <- means$plantcode == code
  subset <- means[index,]
  for (i in 1:dim(subset)[1]){
    new_name <- mC[mC[,1]==code,2]
    #print(!subset[i,]$done)
    if (!subset[i,]$done){
      subset[i,c(5)] <- new_name
      subset[i,c(9)] <- TRUE
      subset$done <- as.logical(subset$done)
    }
  }
  means[index, ] <-  subset
}

means[rownames(means) == "WG4P1FICUCO_P.JPG",]$plantcode <- "FICUCO"
means[rownames(means) == "WG4P2FICUCO_P.JPG",]$plantcode <- "FICUCO"

# The list of codes from the means dataset is shorter than the one from wanang
# specCodesWanang[! (specCodesWanang  %in% as.matrix(sort(unique(means$plantcode))))]

# Strategy
# Run through each line in the 'means' dataset and check if
# it can be fit into datasetW

rows <- dim(means)[1]
added <- rep(FALSE, rows) # Boolean vector for matched entries

for (rownum in 1:rows){
  # Extract the code and plant name from the "means" dataset.
  row <- means[rownum,]  # Subset a single row
  code <- row$code
  #print(code) # From a subset extract the code for a plot e.g. "OG2P3"
  plant <- row$plantcode
  # If the subset is not empty, i.e. no plants with a given name in that plot
  if (length(datasetW[datasetW$CODE == code & datasetW$SP_CODE == plant, ]$AREA) != 0){
    datasetW[datasetW$CODE == code & datasetW$SP_CODE == plant, ]$AREA <- row$meanAreaReal
    datasetW[datasetW$CODE == code & datasetW$SP_CODE == plant, ]$HERB <- row$meanHerb
    added[rownum] <- TRUE
  }
}


###############################################################################
###############################################################################
############################# REPEAT FOR OHU ##################################
######################## BUT WITHOUT CORRECTIONS ##############################
###############################################################################
###############################################################################

# Using only leaf for wanang
# I can repeat the whole procedure for Ohu as well
leaf <- leafOhu
meanHerb <- tapply(leaf$Herbivory, leaf$Label, mean, na.rm = TRUE) 
meanAreaFull <- tapply(leaf$Area.x, leaf$Label, sum, na.rm = TRUE) 
meanAreaReal <- tapply(leaf$Area.y, leaf$Label, sum, na.rm = TRUE) 
means <- as.data.frame(cbind(meanHerb, meanAreaFull, meanAreaReal))
means$filename <- rownames(means)
means$filename <- substr(means$filename, 1, nchar(means$filename)-6) # Drop the _P.jpg part
means$plantcode <- substr(means$filename, 6, nchar(means$filename)) # Extract plant names from photos
means$garden <- substr(means$filename,1, 3) 
means$plot <- substr(means$filename, 4, 5)
means$code <- paste(means$garden, means$plot, sep="")

# First clean the data and then merge them
datasetO <- dataset[dataset$LOCATION == "O",]
# # From the original dataset, what are the species names
# specNamesOhu <- as.matrix(sort(unique(datasetO$SPEC)))
# specCodesOhu <- as.matrix(sort(unique(datasetO$SP_CODE)))
# spComp <- cbind(specCodesOhu,rbind(specNamesOhu,1))

datasetO$SP_CODE <- as.character(datasetO$SP_CODE)

# Calculate the LDMC
datasetO$LDMC <- datasetO$DRY/datasetO$WET..g.

# Get the names from the leaf frame data
#as.matrix(sort(unique(means$plantcode)))
means$plantcode <- as.character(means$plantcode)


# Change the plantcode in means dataset according to "meansCorrect.txt"
# mC <- read.table("meansCorrect.csv", header = FALSE, sep="\t")
# mC[,1] <- as.character(mC[,1])
# mC[,2] <- as.character(mC[,2])

### For loop to correct the codes in 'means' dataset
# means$done <- FALSE
# 
# for (code in mC[,1]){
#   #print(code)
#   index <- means$plantcode == code
#   subset <- means[index,]
#   for (i in 1:dim(subset)[1]){
#     new_name <- mC[mC[,1]==code,2]
#     #print(!subset[i,]$done)
#     if (!subset[i,]$done){
#       subset[i,c(5)] <- new_name
#       subset[i,c(9)] <- TRUE
#       subset$done <- as.logical(subset$done)
#     }
#   }
#   means[index, ] <-  subset
# }

# The list of codes from the means dataset is shorter than the one from wanang
rows <- dim(means)[1]
added <- rep(FALSE, rows) # Boolean vector for matched entries

for (rownum in 1:rows){
  # Extract the code and plant name from the "means" dataset.
  row <- means[rownum,]  # Subset a single row
  code <- row$code
  #print(code) # From a subset extract the code for a plot e.g. "OG2P3"
  plant <- row$plantcode
  #print(plant)# From subset extract plant code as well
  # If the subset is not empty, i.e. no plants with a given name in that plot
  #
  if (length(datasetO[datasetO$CODE == code & datasetO$SP_CODE == plant, ]$AREA) != 0){
    datasetO[datasetO$CODE == code & datasetO$SP_CODE == plant, ]$AREA <- row$meanAreaReal
    datasetO[datasetO$CODE == code & datasetO$SP_CODE == plant, ]$HERB <- row$meanHerb
    added[rownum] <- TRUE
  }
}

# means[!added,]
# datasetO[is.na(datasetO$AREA),]

# Specify which datasets to save:
save <- c("datasetW", "datasetO")
rm(list=ls()[!(ls()%in% save)])


