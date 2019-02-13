# 0. Load data and packages ----

# library(lme4)
# library(lmerTest)
library(vegan)
library(betareg)
library(stats4)
library(ggplot2)

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
main$TREAT <- as.character(main$TREAT)
# subset tree dataset
tree <- main %>%
  filter(LIFE.FORM %in% c("tree", "shrub"))

# Descriptive statistics
test <- read.table("datasets/wng_main_test.txt")
tree_test <- read.table("datasets/wng_tree_test.txt")

# 0.b Remove insecticide ----
test <- test[rownames(test) != "WG1P6", ]
tree_test <- tree_test[rownames(tree_test) != "WG1P6", ]

# Calculate BC dissimilarity 

met <- c("jaccard")
# For each treatment calculate the intra-treatment dissimilarity
treat_str <- as.character(unique(main$TREAT))
disdf <- data.frame()
for(trt in 1:length(treat_str)){
  treat <- treat_str[trt]
  subtreat <- main[main$TREAT == treat, ]
  ctsub <- contingencyTable2(subtreat, "CODE", "SP_CODE", "WEIGHT")
  subdf <- data.frame(bcdis = c(vegdist(ctsub, method=met)), treat = treat)
  disdf <- rbind(disdf, subdf)
}

cibeta <- function(x){
  x <- x[x>0]
  nloglikbeta = function(mu, sig) {
    alpha = mu^2*(1-mu)/sig^2-mu
    beta = alpha*(1/mu-1)
    -sum(dbeta(x, alpha, beta, log=TRUE))
  }
  est = mle(nloglikbeta, start=list(mu=mean(x), sig=sd(x)))
  cest <- confint(est)
  res <- cest[1,]
  mi <- res[1]
  ma <- res[2]
  mu <- est@coef[1]
  return(data.frame(ymin=mi,ymax=ma,y=mu))
}

disdf$treat <- relevel(disdf$treat, "CONTROL")

# png("figs/intradiss.png")
ggplot(disdf,  aes(x = treat, y = bcdis))+
  geom_jitter(width=0.1, color="grey80") + 
  theme_bw() + 
  stat_summary(fun.data=cibeta, color = "grey30",
               geom="pointrange",cex=0.8,lwd=2,
               position = position_dodge(width = 0.90))
# dev.off()
bintra0 <- betareg(bcdis~1, data=disdf)
bintra1 <- betareg(bcdis~treat, data=disdf)

summary(bintra)

# Between treatment dissimilarity
combs <- expand.grid(treat_str,treat_str)
combs <- combs[(combs$Var1 != combs$Var2), ]     
combs <- combs[combs$Var1 == "CONTROL",]

interdf <- data.frame()
for(comb in 1:dim(combs)[1]){
  cb <- c(as.character(combs[comb, ]$Var1), 
          as.character(combs[comb, ]$Var2))
  subdat <- main[main$TREAT %in% cb, ]
  sub_ct <- contingencyTable2(subdat, "CODE", "SP_CODE", "WEIGHT")
  # Treatment data
  plot_tr <- data.frame(tapply(subdat$TREAT, subdat$CODE, unique))
  plot_tr$code <- rownames(plot_tr)
  names(plot_tr) <- c("treat", "code")
  plot_tr <- plot_tr[complete.cases(plot_tr),]
  
  # dissimilarity
  vdbc <- as.matrix(vegdist(sub_ct, method = met))
  inter_comp <- vdbc[plot_tr[plot_tr$treat == cb[1], ]$code,
       plot_tr[plot_tr$treat == cb[2], ]$code]
  subdf <- data.frame(bcvals = c(inter_comp), 
                      cross = paste(cb[1],"-", cb[2], sep = ""))
  interdf <- rbind(interdf, subdf)
}

# png("figs/ctdiss.png")
ggplot(interdf,  aes(x = cross, y = bcvals))+
  geom_jitter(width=0.1, color="grey80") + 
  theme_bw() + 
  stat_summary(fun.data=cibeta, color = "grey30",
               geom="pointrange",cex=0.8,lwd=2,
               position = position_dodge(width = 0.90))  

# dev.off()
bcr <- betareg(bcvals~cross, data=interdf)
summary(bcr)


# library(multcompView)
library(emmeans)

marginal = lsmeans(bcr,~ cross)
CLD = cld(marginal,alpha=0.1,Letters=letters,adjust="tukey")

marginal = lsmeans(bintra1,~ treat)
CLD = cld(marginal,alpha=0.1,Letters=letters,adjust="tukey")

# Balsega beta
library(betapart)
require(vegan)

partition <- function(treat = "CONTROL"){
  con <- main[main$TREAT == treat, ]
  con$SP_CODE <- as.character(con$SP_CODE) 
  test <- contingencyTable2(con, "CODE", "SP_CODE", "WEIGHT")
  return(beta.multi.abund(test, index.family="bray"))
}


pp <- partition()
partition("INSECTICIDE")
partition("FUNGICIDE")
partition("PREDATOR")
partition("WEEVIL125")
partition("WEEVIL25")

beta_part <- data.frame()
for (treat in unique(main$TREAT)){
  res <- partition(treat)
  sub_df2 <- data.frame(treat = treat,
                       type = "Balanced", 
                       val = res$beta.BRAY.BAL)
  sub_df1 <- data.frame(treat = treat,
                       type = "Gradient", 
                       val = res$beta.BRAY.GRA)
  beta_part <- rbind(beta_part, sub_df1, sub_df2)
}

beta_part$treat <- factor(beta_part$treat, levels = c("CONTROL",
                                                      "FUNGICIDE",
                                                      "INSECTICIDE",
                                                      "PREDATOR",
                                                      "WEEVIL25",
                                                      "WEEVIL125"))

# Stacked
# png("figs/beta_part.png", width=400, height=400)
ggplot(beta_part, aes(fill=type, y=val, x=treat)) + 
  geom_bar( stat="identity") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()
