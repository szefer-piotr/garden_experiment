write.table(dataRC_facets, "datasets/dataRC")

read.table("datasets/dataRC", header=T)
getwd()

rcdat <- dataRC_facets[49:60,] # insecticide

# CWM for trees
names(tree)
tree$codetrt <- 
tree$props <- tapply(tree$WEIGHT,)
tree$slaw <- tree$WEIGHT * tree$



sladat <- cwmpanel[cwmpanel$type == "Specific Leaf Area" & cwmpanel$treat %in% c("C","I"), ]
sladat <- cwmpanel[cwmpanel$type == "Leaf Dry Matter Content" & cwmpanel$treat %in% c("C","I"), ]

rcdat$plot <- as.character(rcdat$plot)
sladat$plot <- rownames(sladat)

rcdat$plot == sladat$plot
rcdat <- rcdat[order(rcdat$plot),]
sladat <- sladat[order(sladat$plot),]

rcdat <- rcdat[rcdat$plot != "WG1P6", ]

cordat <- cbind(rcdat, sladat)
cordat[, c(2,6)]

x11(4,4)
plot(meanRC~val, data=cordat)
summary(lm(meanRC~val, data=cordat))
