# Rank abundance curves
main
codes <- as.character(unique(main$CODE))
treats <- c()
for(i in codes){
  treats <- c(treats, as.character(main[main$CODE == i, ]$TREAT[1]))
}

treats <- as.data.frame(treats)
rownames(treats) <- codes
treats$V1 <- as.numeric(treats[,1])

for(i in codes){
  data <- main[main$CODE == i, c("SP_CODE", "WEIGHT")]
  data$sW <- data$WEIGHT/sum(data$WEIGHT)
  plot(data[order(data$sW, decreasing = T), "sW"], type = "l",
       ylim = c(0,1), xlim = c(1,25), col=treats[i,1])
  par(new=T)
}
par(new=F)

library(codyn)
head(collins08)
gdatb1 <- main[main$BLOCK == "WG1", c("TREAT","PLOT", "SP_CODE", "WEIGHT")]
RAC_difference(gdat, species.var = "SP_CODE", 
               abundance_var = "WEIGHT",
               treatment.var = "TREAT", 
               block.var = )



data(pplots)
# With block and no time
df <- subset(pplots, year == 2002 & block < 3)
RAC_difference(df = df,species.var = "species",
               abundance.var = "relative_cover",
               treatment.var ='treatment',
               block.var = "block",
               replicate.var = "plot")
