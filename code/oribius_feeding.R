feed <- read.csv("datasets/oribius.csv", header = T)
names(feed) <- c("species","feeding","genus","family")
feed$species <- as.character(feed$species)
strsplit(feed$species, " ") -> ss

fams <- c()
for(i in 1:length(ss)){fams <- c(fams, ss[[i]][1])}

feed$genus <- fams
specs <- as.vector(as.character(gsub(",","", feed$species)))
feed$species <- specs
specs <- as.vector(as.character(gsub(":","", feed$species)))
feed$species <- specs


library(ggplot2)
library(forcats)
library(RColorBrewer)

# Classic palette BuPu, with 4 colors
coul = brewer.pal(9, "Set1") 

# I can add more tones to this palette :
coul = colorRampPalette(coul)(18)

# Plot it
pie(rep(1, length(coul)), col = coul , main="") 


feed <- feed[order(feed$feeding, decreasing = T),]
p <- ggplot(feed,aes(x=fct_reorder(species, feeding,.desc =TRUE), 
                      y=feeding, fill = family))

png("figs/figS1.png", height=10, width=10,units = "in", res=400)
p + geom_bar(stat = "identity") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        axis.title.x=element_blank()) + ylab("Feeding") +
  scale_fill_manual(values = coul, name="Family")
dev.off()


