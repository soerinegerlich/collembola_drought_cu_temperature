# looking at growth as the response variable

# Load library
library(tidyverse) # v. 2.0.0
library(patchwork)
library(emmeans)


# Input origin data
dat.org <- read_csv('data/sugrwe.csv') %>% as_tibble() %>% mutate(., vict = number-survivor)

# check what the data looks like

summary(dat.org$growth)

# problem with right skewness?

sum(dat.org$growth == 0, na.rm = TRUE)

mean(dat.org$growth == 0, na.rm = TRUE) * 100

hist(dat.org$growth,
     breaks = 30,
     main = "Body weight distribution",
     xlab = "Body weight (mg)")

boxplot(dat.org$growth,
        ylab = "Body weight (mg)")

dat.org %>%
  summarise(
    n = n(),
    zeros = sum(growth == 0, na.rm = TRUE),
    pct_zeros = 100 * mean(growth == 0, na.rm = TRUE),
    min = min(growth, na.rm = TRUE),
    q1 = quantile(growth, 0.25, na.rm = TRUE),
    median = median(growth, na.rm = TRUE),
    mean = mean(growth, na.rm = TRUE),
    q3 = quantile(growth, 0.75, na.rm = TRUE),
    max = max(growth, na.rm = TRUE)
  )

# does log() help

hist(log(dat.org$growth))
qqnorm(log(dat.org$growth))
qqline(log(dat.org$growth))

mod.growth <- lm(
  log(growth) ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.org
)

par(mfrow = c(2,2))
plot(mod.growth)

mod.growth.glm <- glm(
  growth ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  family = Gamma(link = "log"),
  data = dat.org
)


par(mfrow = c(2,2))
plot(mod.growth.glm)

###################### separate in cold and warm ##############################

dat.cold <- dat.org %>%
  filter(temp %in% c(20.1, 9.7, 4.5))


dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1, 27.7))

#dat.warm <- dat.warm[-c(161, 169, 237), ]

# cold

mod.growth.cold <- lm(
  log(growth) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  data = dat.cold
)

par(mfrow = c(2,2))
plot(mod.growth.cold)

# there's an issue with some outliers, but the conclusion of the analysis does
# not change when removing the outliers so they are kept

mod.growth.glm.cold <- glm(
  growth ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = Gamma(link = "log"),
  data = dat.cold
)

summary(mod.growth.glm.cold)

par(mfrow = c(2,2))
plot(mod.growth.glm.cold)

hist(dat.cold$growth,
     breaks = 30,
     main = "Body weight distribution",
     xlab = "Body weight (mg)")

hist(sqrt(dat.cold$growth),
     breaks = 30,
     main = "Body weight distribution",
     xlab = "Body weight (mg)")

# square root of the response variable does not help

# warm

mod.growth.warm <- lm(
  log(growth) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  data = dat.warm
)

par(mfrow = c(2,2))
plot(mod.growth.warm)

mod.growth.glm.warm <- glm(
  growth ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = Gamma(link = "log"),
  data = dat.warm
)
summary(mod.growth.glm.warm)

par(mfrow = c(2,2))
plot(mod.growth.glm.warm)

hist(dat.warm$growth,
     breaks = 30,
     main = "Body weight distribution",
     xlab = "Body weight (mg)")

hist(log(dat.warm$growth),
     breaks = 30,
     main = "Body weight distribution",
     xlab = "Body weight (mg)")

# as with the cold gradient, there's an issue with some outliers, 
# but the conclusion of the analysis does not change when removing the 
# outliers so they are kept

# lm with log transformed growth response variable looks as good as the glm.
# the choice is then to use lm() because it is a simpler model and the output
# is easier to understand

mod.growth.cold.full <- lm(
  log(growth) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  data = dat.cold
)

summary(mod.growth.cold.full)

drop1(mod.growth.cold.full, test = "F")

library(effectsize)

eta_squared(mod.growth.cold.full, partial = TRUE)

library(car)

anova_growth <- Anova(mod.growth.cold.full, type = 3)

effectsize::eta_squared(
  anova_growth,
  partial = TRUE,
  ci = NULL
)

mod.growth.cold.twoway <- lm(
  log(growth) ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.cold
)

summary(mod.growth.cold.twoway)

mod.growth.cold.simple <- lm(
  log(growth) ~ factor(temp) + factor(lgswp) + factor(sqrcu),
  data = dat.cold
)

summary(mod.growth.cold.simple)

anova(mod.growth.cold.twoway, mod.growth.cold.simple)

drop1(mod.growth.cold.twoway, test = "F")

mod.growth.cold.final <- lm(
  log(growth) ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp) +
    factor(temp):factor(sqrcu),
  data = dat.cold
)

drop1(mod.growth.cold.final, test = "F")

par(mfrow = c(2,2))
plot(mod.growth.cold.final)

anova(mod.growth.final, mod.growth.cold.full)

drop1(mod.growth.final, test = "F")

# the final model includes temperature interactions with 
# SWP and Cu contamination.

# warm

mod.growth.warm.full <- lm(
  log(growth) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  data = dat.warm
)

summary(mod.growth.warm.full)

drop1(mod.growth.warm.full, test = "F")

library(effectsize)

eta_squared(mod.growth.warm.full, partial = TRUE)

library(car)

anova_growth <- Anova(mod.growth.warm.full, type = 3)

effectsize::eta_squared(
  anova_growth,
  partial = TRUE,
  ci = NULL
)

mod.growth.warm.twoway <- lm(
  log(growth) ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.warm
)


summary(mod.growth.warm.twoway)

drop1(mod.growth.warm.twoway, test = "F")

mod.growth.warm.simple <- lm(
  log(growth) ~ factor(temp) + factor(lgswp) + factor(sqrcu),
  data = dat.warm
)

summary(mod.growth.warm.simple)

anova(mod.growth.warm.full, mod.growth.warm.twoway, mod.growth.warm.simple)

drop1(mod.growth.warm.twoway, test = "F")


par(mfrow = c(2,2))
plot(mod.growth.warm.twoway)


# the model containing all two way interactions is the best: mod.growth.warm.twoway

# see the script: growth_figures for creating figure 4 in the main manuscript


##########


