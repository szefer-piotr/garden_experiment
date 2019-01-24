#Dataset preparation
dataset <- read.csv("datasets/supplementary/main_biomass_wng_ohu.csv")
dataset$LOCATION <- substr(dataset$CODE, 1, 1)

#Subsetting for Wanang
dataset <- dataset[dataset$LOCATION == "W", ]
dataset$CODE <- factor(dataset$CODE)
dataset$WEIGHT <- as.numeric(as.character(dataset$WEIGHT))

# Treat is a subset with GARDEN, PLOT, PLOT_CODE, TREATMENT VALUES.
treat <- data.frame(PLOT_CODE = NULL, TREATMENT=NULL)
for (plotCode in unique(dataset$CODE)){
  Row <- data.frame(PLOT_CODE = plotCode, TREATMENT = dataset[dataset$CODE == plotCode, ]$TREAT[1])
  treat <- rbind(treat, Row)
}
treat$GARDEN <- substr(as.character(treat$PLOT_CODE), 1,3)
treat$PLOT <- substr(as.character(treat$PLOT_CODE), 4,5)

plants <- table(dataset$CODE, dataset$SP_CODE)
plants <- plants[rowSums(plants) != 0, colSums(plants) != 0]
allSpecCodes <- colnames(plants)
allPlotCodes <- rownames(plants)
entries <- which(plants != 0, arr.ind = TRUE)


for (entry in 1:dim(entries)[1]){
  plot <- entries[entry,1]
  plant <- entries[entry,2]
  
  specCode <- allSpecCodes[plant]
  plotCode <- allPlotCodes[plot]

  plants[plot,plant] <- sum(dataset[dataset$CODE == plotCode & dataset$SP_CODE == specCode,]$WEIGHT)
}
plants[is.na(plants)] <- 0

# Treatment dataset
treatment2 <- treat
treatment2$CON <- ifelse(treat$TREAT == "CONTROL", 1, 0)
treatment2$FUN <- ifelse(treat$TREAT == "FUNGICIDE", 1, 0)
treatment2$HLO <- ifelse(treat$TREAT == "WEEVIL25", 1, 0)
treatment2$HHI <- ifelse(treat$TREAT == "WEEVIL125", 1, 0)
treatment2$INS <- ifelse(treat$TREAT == "INSECTICIDE", 1, 0)
treatment2$PRE <- ifelse(treat$TREAT == "PREDATOR", 1, 0)
treatment2 <- treatment2[order(as.character(treatment2$PLOT_CODE)), ]
# Analysis


# Full Analysis
library(vegan)
RDA1d <- rda(log(plants+1)~FUN+HLO+HHI+INS+PRE+Condition(GARDEN),
             data = treatment2)
summary(RDA1d)
anova(RDA1d, by="terms", permutations = 9999)

# Plotting

# Getting colours ready
cols <- c(rgb(155,  0,  0,255,maxColorValue = 255),
          rgb(  0,155,  0,255,maxColorValue = 255),
          rgb(  0,  0,155,255,maxColorValue = 255),
          rgb(  0,155,155,255,maxColorValue = 255),
          rgb(155,155,  0,255,maxColorValue = 255),
          rgb(155,  0,155,255,maxColorValue = 255))
numTreat <- as.numeric(treatment2$TREAT)
pointChars = as.numeric(treatment2$TREAT)
colors <- numTreat
for (col in 1:6){
  colors[colors == col] <- cols[col]
}

# Main plot
plot(RDA1d, display = "sites", type = "n", 
     xlab = "RDA1 [12.56%]", ylab="RDA2 [4.41%]") 
#xlim = c(-5, 5), ylim = c(-5,5))
points(RDA1d, pch = pointChars , cex = 1.5, col = colors, lwd = 2)
text(RDA1d, display="bp", lwd = 2,col = c(cols[2], cols[6],
                                          cols[5], cols[3],
                                          cols[4]))

legend('bottomright', c("CONTROL", "FUNGICIDE", "HERBIVORY lev. 1",
                        "HERBIVORY lev. 2", "INSECTICIDE", "PREDATOR"), 
       cex = 0.7,
       pch = c(1, 2, 6, 5, 3, 4), 
       col = c(cols[1],cols[2], cols[6],cols[5], cols[3],cols[4]))


##### GOODNESS OF FIT ANALYSIS
plot(RDA1d, display = "species", type = "n", scaling = 2.5) 
text(RDA1d, display="bp", lwd = 2,col = c(cols[2], cols[6],
                                          cols[5], cols[3],
                                          cols[4]))
scores <- as.data.frame(scores(RDA1d, display="sp"))
TRESHOLD <- 0.05

GOOD <- goodness(RDA1d)[,c(1,2)]
GOOD <- GOOD[GOOD[,1] > TRESHOLD | GOOD[,1] > TRESHOLD, ]
scores <- scores[rownames(GOOD),]

scaling <- 1
labels <- rownames(scores)
gap <- 1.1
for (species in 1:length(labels)){
  endpoints <- as.vector(scores[labels[species],])
  endpoints <- as.numeric(endpoints)*scaling
  #lines(c(0,endpoints[1]), c(0,endpoints[2]), col = "lightgray")
  arrows(0,0,endpoints[1],endpoints[2], col = "lightgray", length = 0.05)
  lab <- as.character(labels[species])
  #text(endpoints[1]*gap, endpoints[2]*gap, lab, cex = 0.5,col = "gray50")
}

selection <- rownames(summary(RDA1d)$species) %in% rownames(GOOD)
orditorp(RDA1d, display = "species", select = selection, col = "red")
orditorp(RDA1d, display = "species")


names(summary(RDA1d))
