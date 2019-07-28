# Cwm analysis

library(dplyr)

# full community ----
main <- read.table("datasets/wng_main_clean.txt", header = T) #506 rows


main <- main %>%
  filter(!is.na(main$SLA) & main$SLA < 1000 & main$SLA != 0) #496 obs

test <- tapply(main$WEIGHT, main$CODE, 
               function(X){return(X/sum(X))}) #497 obs

test <- stack(test)
test$sla <- main$SLA # here!
test$ldmc <- main$LDMC
test$herb <- main$HERB

test$garden <- substr(test$ind, 1,3)
test$plot <- substr(test$ind, 4,5)
test$treat <- 0

for (row in 1:dim(test)[1]) {
  code <- test$ind[row]
  to_assign <- as.character(main[main$CODE == code,]$TREAT[1])
  test[row,]$treat <- to_assign 
}

# get cwm sla values
names(test) <- c("prop", "code", "sla", "ldmc_logit", "herb_logit", "garden","plot","treat")
head(test)
test$propsla <- test$prop * test$sla
test$propldmc <- test$prop * test$ldmc_logit
test$propherb <- test$prop * test$herb_logit


# group variables
data <- data.frame(cwmsla = tapply(test$propsla, test$code, sum))
data$cwmldmc = tapply(test$propldmc, test$code, sum)
data$cwmherb = tapply(test$propherb, test$code, sum)

data$code <- rownames(data)
data$garden <- substr(data$code, 1,3)
data$plot <- substr(data$code, 4,5)
test$treat <- as.character(test$treat)
data$treat <- tapply(test$treat, test$code, unique)
traitdat <- data[data$code != "WG1P6", ]

# tests
tapply(test$prop, test$code, sum)

# >>>>Stat Tests -----
library(ggplot2)
p <- ggplot(traitdat, aes(x = treat, y = cwmsla))
p + geom_point()

c <- ggplot(traitdat, aes(x = treat, y = cwmldmc))
c + geom_point()

h <- ggplot(traitdat, aes(x = treat, y = cwmherb))
h + geom_point()

library(lme4)
library(lmerTest)
sla_rbl <- lmer(cwmsla ~ treat + (1|garden), data = traitdat)
ldmc_rbl <- lmer(cwmldmc ~ treat + (1|garden), data = traitdat)
herb_rbl <- lmer(cwmherb ~ treat + (1|garden), data = traitdat)

summary(sla_rbl)
summary(ldmc_rbl)
summary(herb_rbl)

# >>>>Plots ------

# Get the data into right format
cwmpanel <- rbind(setNames(traitdat[,c("treat", "cwmsla" )], c("treat", "val")),
                  setNames(traitdat[,c("treat", "cwmldmc")], c("treat", "val")))
cwmpanel$type <- rep(c("Specific Leaf Area", "Leaf Dry Matter Content"), 
                     each = dim(traitdat)[1])
cwmpanel$treat <- factor(cwmpanel$treat, labels = c("F","H2","C","P","H1","I"))
cwmpanel$treat <- ordered(cwmpanel$treat, labels = c("C","F","I","P","H1","H2"))

  
fig2 <- ggplot(cwmpanel, aes(x=treat, y=val)) + 
  facet_wrap(~type, scales = "free") + 
  geom_point(col = "grey80", alpha = 0.5, cex = 4) + 
  theme_bw()


colors = color=c("black","grey50","red","grey50","grey50","grey50",
                 "black","grey50","orange","grey50","grey50","grey50")

# png("figs/fig2.png",width=800, height = 400)
png("figs/fig2.png",width=8, height = 4,units = 'in',res=1200)

# fig2 + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
#                     geom="errorbar", color= colors, width=0.2, lwd=1.5) +
#   stat_summary(fun.y=mean, geom="point", color=colors, cex = 5) +
#   theme(axis.text.x=element_text(angle=0, size=20, hjust=0.5),
#         axis.text.y=element_text(angle=0, size=20, hjust=0.5),
#         strip.text = element_text(size=20),
#         legend.justification=c(0.5,0.5), 
#         legend.position="bottom")+xlab("")+ylab("")

fig2 + stat_summary(fun.data=mean_cl_boot, 
                    geom="pointrange", color= colors, width=0.2, lwd=1.5) +
  stat_summary(fun.y=mean, geom="point", color=colors, cex = 5) +
  theme(axis.text.x=element_text(angle=0, size=20, hjust=0.5),
        axis.text.y=element_text(angle=0, size=20, hjust=0.5),
        strip.text = element_text(size=20),
        legend.justification=c(0.5,0.5), 
        legend.position="bottom")+xlab("")+ylab("")


dev.off()

# tree community ----

main <- read.table("datasets/wng_main_clean.txt", header = T) #506 rows


main <- main %>%
  filter(!is.na(main$SLA) & main$SLA < 1000 & main$SLA != 0) %>%#496 obs
  filter(LIFE.FORM %in% c("tree", "shrub"))

test <- tapply(main$WEIGHT, main$CODE, 
               function(X){return(X/sum(X))}) #497 obs

test <- stack(test)
test$sla <- main$SLA # tutaj!
test$ldmc <- main$LDMC
test$herb <- main$HERB

test$garden <- substr(test$ind, 1,3)
test$plot <- substr(test$ind, 4,5)
test$treat <- 0

for (row in 1:dim(test)[1]) {
  code <- test$ind[row]
  to_assign <- as.character(main[main$CODE == code,]$TREAT[1])
  test[row,]$treat <- to_assign 
}

# get cwm sla values
names(test) <- c("prop", "code", "sla", "ldmc_logit", "herb_logit", "garden","plot","treat")
head(test)
test$propsla <- test$prop * test$sla
test$propldmc <- test$prop * test$ldmc_logit
test$propherb <- test$prop * test$herb_logit


# group variables
data <- data.frame(cwmsla = tapply(test$propsla, test$code, sum))
data$cwmldmc = tapply(test$propldmc, test$code, sum)
data$cwmherb = tapply(test$propherb, test$code, sum)

data$code <- rownames(data)
data$garden <- substr(data$code, 1,3)
data$plot <- substr(data$code, 4,5)
test$treat <- as.character(test$treat)
data$treat <- tapply(test$treat, test$code, unique)
traitdat <- data[data$code != "WG1P6", ]

# tests
tapply(test$prop, test$code, sum)

# >>>>Stat Tests -----
library(ggplot2)
library(lmerTest)
p <- ggplot(data, aes(x = treat, y = cwmsla))
p + geom_point()

c <- ggplot(data, aes(x = treat, y = cwmldmc))
c + geom_point()

h <- ggplot(data, aes(x = treat, y = cwmherb))
h + geom_point()

library(lme4)
sla_rbl <- lmer(cwmsla ~ treat + (1|garden), data = traitdat)
ldmc_rbl <- lmer(cwmldmc ~ treat + (1|garden), data = traitdat)
herb_rbl <- lmer(cwmherb ~ treat + (1|garden), data = traitdat)

summary(sla_rbl)
summary(ldmc_rbl)
summary(herb_rbl) # marginal predator and weevil
