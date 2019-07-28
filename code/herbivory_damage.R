# Herbivory damage plot
library(lme4)
library(ggplot2)
library(lmerTest)

main <- read.table("datasets/wng_main_clean.txt", header = T)

# Remoove NA's
main <- main[!is.na(main$HERB), ]

# Remove zero's
mainz <- main[main$HERB == 0, ]
main <- main[main$HERB < 0, ]

# mgras <- main[main$LIFE.FORM %in% c("grass"), ]
# msedg <- main[main$LIFE.FORM %in% c("sedge"), ]
# mvine <- main[main$LIFE.FORM %in% c("vine"), ]
# mherb <- main[main$LIFE.FORM %in% c("herb"), ]
# mtree <- main[main$LIFE.FORM %in% c("tree"), ]
# mshrub <- main[main$LIFE.FORM %in% c("shrub"), ]
# mwood <- main[main$LIFE.FORM %in% c("tree", "shrub"), ]
# 
# # Plot the data
# plt <- ggplot(mgras, aes(x = TREAT, y = HERB))
# plt + geom_jitter(width =0.1)
# 
# plt <- ggplot(msedg, aes(x = TREAT, y = HERB))
# plt + geom_jitter(width =0.1)
# 
# plt <- ggplot(mvine, aes(x = TREAT, y = HERB))
# plt + geom_jitter(width =0.1)
# 
# plt <- ggplot(mherb, aes(x = TREAT, y = HERB))
# plt + geom_jitter(width =0.1)
# 
# plt <- ggplot(mwood, aes(x = TREAT, y = HERB))
# plt + geom_jitter(width =0.1)
# 
# # Main with types
# grasmod <- lm(HERB ~ TREAT, data = mgras)
# summary(grasmod)
# sedgmod <- lm(HERB ~ TREAT, data = msedg)
# summary(sedgmod)
# vinemod <- lm(HERB ~ TREAT, data = mvine)
# summary(vinemod)
# herbmod <- lm(HERB ~ TREAT, data = mherb)
# summary(herbmod)
# woodmod <- lm(HERB ~ TREAT, data = mwood)
# summary(woodmod)
# treemod <- lm(HERB ~ TREAT, data = mtree)
# summary(treemod)
# shrubmod <- lm(HERB ~ TREAT, data = mshrub)
# summary(shrubmod)

# Simple LM - no random effects
mod1 <- lm(HERB ~ TREAT, data = main)
summary(mod1)

hpl <- ggplot(main, aes(x = TREAT, y = HERB))
hpl + geom_jitter(aes(colour = LIFE.FORM), width =0.2, alpha= 0.5) + 
  xlab("") + 
  ylab("logit[Prop. of area missing]") +
  stat_summary(fun.data=mean_cl_boot, 
               geom="pointrange", 
               col = "grey30")

# hpl + geom_jitter(aes(colour = LIFE.FORM), width =0.2, alpha= 0.5) + 
#   xlab("") + 
#   ylab("logit[Prop. of area missing]") +
#   stat_summary(fun.data=smean_cl_normal, 
#                geom="errorbar", 
#                col = "grey30",
#                width=0.2, lwd=1) + 
#   stat_summary(geom= "point", fun.y=mean)


