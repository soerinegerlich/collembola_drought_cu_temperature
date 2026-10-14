# reproduction figure for main manuscript

# looking at reproduction as the response variable

# Load library
library(tidyverse) # v. 2.0.0
library(patchwork)
library(emmeans)
library(MASS)

# Input origin data
dat.org <- read_csv('data/sugrwe.csv') %>% as_tibble() %>% mutate(., vict = number-survivor)

# check data distribution and issue with zeros

summary(dat.org$repro)

boxplot(repro ~ temp, data = dat.org)

aggregate(repro ~ temp,
          data = dat.org,
          mean)

with(dat.org,
     tapply(repro,
            temp,
            function(x) mean(x == 0)))

# first check glm with a poisson distribution

mod.pois <- glm(
  repro ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = poisson,
  data = dat.org
)

deviance(mod.pois) /
  df.residual(mod.pois)

# much larger than 1 and Poisson is therefor inappropriate

# check negative binomial

mod.nb <- glm.nb(
  repro ~ factor(temp) *
    factor(lgswp) *
    factor(sqrcu),
  data = dat.org
)

plot(mod.nb)

################################################################################
# Only look at 20.1 and 24.1 degree C data #
################################################################################

dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1))

# first check glm with a poisson distribution

mod.pois <- glm(
  repro ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = poisson,
  data = dat.warm
)

deviance(mod.pois) /
  df.residual(mod.pois)

# much larger than 1 and Poisson is therefor inappropriate

# check negative binomial

mod.nb <- glm.nb(
  repro ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  data = dat.warm
)

summary(mod.nb)

drop1(mod.nb, test = "F")

mod.offspring.final <- glm.nb(
  repro ~ (factor(temp) +
    factor(lgswp) +
    factor(sqrcu))^2,
  data = dat.warm)

summary(mod.offspring.final)

drop1(mod.offspring.final, test = "F")


# remove swp and cu interaction

mod.offspring.final <- glm.nb(
  repro ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp) +
    factor(temp):factor(sqrcu),
  data = dat.warm
)

summary(mod.offspring.final)

drop1(mod.offspring.final, test = "F")

library(DHARMa)

# Simulate residuals
sim.repro <- simulateResiduals(
  fittedModel = mod.offspring.final,
  n = 1000
)

# Standard diagnostic plot
plot(sim.repro)

# Formal tests
testDispersion(sim.repro)
testZeroInflation(sim.repro)
testUniformity(sim.repro)

# Outlier test
testOutliers(sim.repro)


