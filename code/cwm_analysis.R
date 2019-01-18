# Unfinished! Is this a good dataset?

main <- read.table("datasets/wng_main_clean.txt", header = T)

main <- main %>%
  filter(!is.na(main$SLA) & main$SLA < 1000 & main$SLA != 0)

test <- tapply(main$WEIGHT, main$CODE, 
               function(X){return(X/sum(X))})

test <- stack(test)
test$sla <- main$SLA
test$ldmc <- main$LDMC

test$garden <- substr(test$ind, 1,3)
test$plot <- substr(test$ind, 4,5)
test$treat <- main$TREAT

# get cwm sla values
names(test) <- c("prop", "code", "sla", "ldmc_logit", "garden","plot","treat")
head(test)
test$propsla <- test$prop * test$sla
test$propldmc <- test$prop * test$ldmc_logit

# group variables
data <- data.frame(cwmsla = tapply(test$propsla, test$code, sum))
data$code <- rownames(data)
data$garden <- substr(data$code, 1,3)
data$plot <- substr(data$code, 4,5)
test$treat <- as.character(test$treat)
data$treat <- tapply(test$treat, test$code, unique)

# tests



p <- ggplot(test, aes(x = treat, y = p))


ggplot(tgc, aes(x=dose, y=len, colour=supp)) + 
  geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd)
