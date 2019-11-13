# Rank abundance curves
main
codes <- as.character(unique(main$CODE))
treats <- c()
for(i in codes){
  treats <- c(treats, as.character(main[main$CODE == i, ]$TREAT[1]))
}

treats <- as.data.frame(treats)
rownames(treats) <- codes
treats$V1 <- as.numeric(treats$V1)

for(i in codes){
  data <- main[main$CODE == i, c("SP_CODE", "WEIGHT")]
  data$sW <- data$WEIGHT/sum(data$WEIGHT)
  plot(data[order(data$sW, decreasing = T), "sW"], type = "l",
       ylim = c(0,1), xlim = c(1,25), col=treats[i,1], lty=treats[i,1])
  par(new=T)
}
par(new=F)
