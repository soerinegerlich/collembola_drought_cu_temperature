
# full factorial design in modelling internal Cu and Pb concentrations

library(tidyverse) # v. 2.0.0
library(robustbase) # v. 0.99-7

# retrieve dataset

cu.list <- c(17, 436, 1629)
swp.list <- c(3.2, 5, 10, 39.8, 63.1, 100)
temp.list <- c(4.5, 9.7, 20.1, 24.1, 27.7)

dat.org <- read_csv('data/cuaspb.csv') %>% as_tibble() %>% 
  mutate(., 
         sqrcu = case_when(
           cul == 'L1' ~ cu.list[1],
           cul == 'L2' ~ cu.list[2],
           cul == 'L3' ~ cu.list[3]
         ),
         tem = case_when(
           temp == '5' ~ temp.list[1],
           temp == '10' ~ temp.list[2],
           temp == '20' ~ temp.list[3],
           temp == '24' ~ temp.list[4],
           temp == '28' ~ temp.list[5]
         ),
         lgswp = case_when(
           swp == '3' ~ swp.list[1],
           swp == '5' ~ swp.list[2],
           swp == '10' ~ swp.list[3],
           swp == '40' ~ swp.list[4],
           swp == '70' ~ swp.list[5],
           swp == '100' ~ swp.list[6]
         )) %>% dplyr::select(tem, sqrcu, lgswp, weight, cu, pb, repli) %>% 
  filter(., between(cu, 0, 2000) & between(pb, 0, 2000))


# Because the data is not normally distributed, transformation is needed. 
# Quick test to check distribution of residuals
summary(dat.org$cu)
hist(dat.org$cu)
boxplot(dat.org$cu)

# the distribution is extremely right skewed.

mod.log <- lm(
  log(cu + 1) ~ factor(tem) * factor(lgswp) * factor(sqrcu),
  data = dat.org
)

par(mfrow = c(2,2))
plot(mod.log)

# check robust regression

library(lmtest)

bptest(mod.log)

mod_rob <- lmrob(
  log(cu + 1) ~ tem * lgswp * sqrcu,
  data = dat.org
)

par(mfrow = c(2,2))
plot(mod_rob)

hist(log(dat.org$cu))
hist(residuals(mod.log))

###############################################################################
# approach: use log1p() and lm(). The diagnostic plots look fine.
###############################################################################

# change predictors to categorical variables

dat.org$temp_f <- factor(dat.org$tem)
dat.org$swp_f <- factor(dat.org$lgswp)
dat.org$cul_f <- factor(dat.org$sqrcu)

mod_cu <- lmrob(
  log1p(cu) ~ tem * lgswp * sqrcu,
  data = dat.org
)

summary(mod_cu)


mod_cu_f <- lmrob(
  log1p(cu) ~ temp_f * swp_f * cul_f,
  data = dat.org
)

summary(mod_cu_f)

# the model expands to:
#  Temp +
#  SWP +
#  CuLevel +
#  Temp:SWP +
#  Temp:CuLevel +
#  SWP:CuLevel +
#  Temp:SWP:CuLevel

# which matches the factorial design

# however, we want to divide the analysis in cold and warm temp gradients

###############################################################################
#################### cold and warm temperature gradients ######################
###############################################################################

# cold gradient

dat.cold <- dat.org %>%
  filter(tem %in% c(20.1, 9.7, 4.5))

# tests with both robust linear regression and linear regression

mod.cold.full.rob <- lmrob(
  log1p(cu) ~ factor(tem) * factor(lgswp) * factor(sqrcu),
  data = dat.cold
)
# model does not converge

mod.cold.twoway <- lmrob(
  log1p(cu) ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.cold
)

summary(mod.cold.twoway)

# testing linear models

mod.cold.full <- lm(
  log1p(cu) ~ factor(tem) * factor(lgswp) * factor(sqrcu),
  data = dat.cold
)

summary(mod.cold.full)

mod.cold.lm <- lm(
  log1p(cu) ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.cold
)

summary(mod.cold.lm)

drop1(mod.cold.full, test = "F")
drop1(mod.cold.lm, test = "F")

par(mfrow = c(2,2))
plot(mod.cold.lm)

# calculate effect size measures to estimate how much variance in the response
# variables is accounted for by the explanatory variables


library(effectsize)

eta_squared(mod.cold.full, partial = TRUE)

library(car)

anova_cu <- Anova(mod.cold.full, type = 3)

effectsize::eta_squared(
  anova_cu,
  partial = TRUE,
  ci = NULL
)
# warm gradient

dat.warm <- dat.org %>%
  filter(tem %in% c(20.1, 24.1, 27.7))

mod.warm.full <- lm(
  log1p(cu) ~ factor(tem) * factor(lgswp) * factor(sqrcu),
  data = dat.cold
)

summary(mod.warm.full)

drop1(mod.warm.full, test = "F")

mod.warm.twoway <- lmrob(
  log1p(cu) ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.warm
)

# look at results from a normal lm or simplify model by not including 
# interactions

summary(mod.warm.twoway)

mod.warm.lm <- lm(
  log1p(cu) ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.warm
)

summary(mod.warm.lm)

drop1(mod.warm.lm, test = "F")

par(mfrow = c(2,2))
plot(mod.warm.lm)

# calculate effect size measures

eta_squared(mod.warm.lm, partial = TRUE)

anova_cu <- Anova(mod.warm.lm, type = 3)

effectsize::eta_squared(
  anova_cu,
  partial = TRUE,
  ci = NULL
)


# see script "cu_pb_figures" for figures

################################################################################
################## check for pb ################################################
################################################################################

# run full model initially

summary(dat.org$pb)
hist(dat.org$pb)
boxplot(dat.org$pb)

mod.pb.full <- lmrob(
  log(pb + 1) ~ factor(tem) * factor(lgswp) * factor(sqrcu),
  data = dat.org
)

summary(mod.pb.full)

mod.log <- lm(
  log(pb + 1) ~ factor(tem) * factor(lgswp) * factor(sqrcu),
  data = dat.org
)

summary(mod.log)

par(mfrow = c(2,2))
plot(mod.pb.full)

par(mfrow = c(2,2))
plot(mod.log)

mod_rob <- lmrob(
  pb ~ tem * lgswp * sqrcu,
  data = dat.org
)

par(mfrow = c(2,2))
plot(mod_rob)

hist(log(dat.org$cu))
hist(residuals(mod.log))

###############################################################################
# approach: use log1p() Keeps the same approach as before.
# The diagnostic plots look fine.
###############################################################################


dat.org$temp_f <- factor(dat.org$tem)
dat.org$swp_f <- factor(dat.org$lgswp)
dat.org$cul_f <- factor(dat.org$sqrcu)

mod_pb <- lmrob(
  log1p(pb) ~ tem * lgswp * sqrcu,
  data = dat.org
)

summary(mod_pb)

# the model cannot converge

mod2_pb <- lmrob(
  log1p(pb) ~ (tem + lgswp + sqrcu)^2,
  data = dat.org
)

summary(mod2_pb)

# try with factors instead
mod2_pb_f <- lmrob(
  log1p(pb) ~ (temp_f + swp_f + cul_f)^2,
  data = dat.org
)

summary(mod2_pb_f)

###############################################################################
#################### cold and warm temperature gradients ######################
###############################################################################

# cold gradient

dat.cold <- dat.org %>%
  filter(tem %in% c(20.1, 9.7, 4.5))

mod.cold.full <- lm(
  log1p(pb) ~ factor(tem) * factor(lgswp) * factor(sqrcu),
  data = dat.cold
)

summary(mod.cold.full)

mod.cold.twoway <- lmrob(
  log1p(pb) ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.cold
)

summary(mod.cold.twoway)

mod.cold.lm <- lm(
  log1p(pb) ~ factor(tem) + factor(lgswp) + factor(sqrcu),
  data = dat.cold
)

drop1(mod.cold.lm, test = "F")
summary(mod.cold.lm)

par(mfrow = c(2,2))
plot(mod.cold.lm)

# calculate effect size measures on full model to determine effect of 
# three way interaction

library(effectsize)

eta_squared(mod.cold.full, partial = TRUE)

library(car)

anova_cu <- Anova(mod.cold.full, type = 3)

effectsize::eta_squared(
  anova_cu,
  partial = TRUE,
  ci = NULL
)


# We can suffice with main effects and no interactions

# warm gradient

dat.warm <- dat.org %>%
  filter(tem %in% c(20.1, 24.1, 27.7))

mod.warm.full <- lm(
  log1p(pb) ~ (factor(tem) * factor(lgswp) * factor(sqrcu)),
  data = dat.warm
)

summary(mod.warm.full)

mod.warm.twoway <- lm(
  log1p(pb) ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.warm
)

summary(mod.warm.lm)

mod.warm.lm <- lm(
  log1p(pb) ~ factor(tem) + factor(lgswp) + factor(sqrcu),
  data = dat.warm
)

drop1(mod.warm.lm, test = "F")

summary(mod.warm.lm)

par(mfrow = c(2,2))
plot(mod.warm.lm)

# calculate effect size measures

library(effectsize)

eta_squared(mod.warm.full, partial = TRUE)

library(car)

anova_cu <- Anova(mod.warm.full, type = 3)

effectsize::eta_squared(
  anova_cu,
  partial = TRUE,
  ci = NULL
)

# see script "cu_pb_figures" for figures
