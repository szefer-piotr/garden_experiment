# Simualte biomass -----

rm(list = ls())
set.seed(2345)

N <- 30
unit.df <- data.frame(unit = c(1:N), a = rnorm(N))

head(unit.df, 3)
unit.df <-  within(unit.df, {
  E.alpha.given.a <-  1 - 0.15 * a
  E.beta.given.a <-  3 + 0.3 * a
})
head(unit.df, 3)

library(mvtnorm)
q = 0.2
r = 0.9
s = 0.5
cov.matrix <- matrix(c(q^2, r * q * s, r * q * s, s^2), nrow = 2,
                     byrow = TRUE)
random.effects <- rmvnorm(N, mean = c(0, 0), sigma = cov.matrix)
unit.df$alpha <- unit.df$E.alpha.given.a + random.effects[, 1]
unit.df$beta <- unit.df$E.beta.given.a + random.effects[, 2]
head(unit.df, 3)

J <- 30
M = J * N  #Total number of observations
x.grid = seq(-4, 4, by = 8/J)[0:30]

within.unit.df <-  data.frame(unit = sort(rep(c(1:N), J)), j = rep(c(1:J),
                                                                   N), x =rep(x.grid, N))
flat.df = merge(unit.df, within.unit.df)

flat.df <-  within(flat.df, y <-  alpha + x * beta + 0.75 * rnorm(n = M))
simple.df <-  flat.df[, c("unit", "a", "x", "y")]
head(simple.df, 3)

library(lme4)
my.lmer <-  lmer(y ~ x + (1 + x | unit), data = simple.df)
cat("AIC =", AIC(my.lmer))
my.lmer <-  lmer(y ~ x + a + x * a + (1 + x | unit), data = simple.df)
summary(my.lmer)
