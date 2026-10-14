
# modelling survival of springtails in a full factorial design

# Load library
library(tidyverse) # v. 2.0.0
library(patchwork)
library(emmeans)


# Input origin data
dat.org <- read_csv('data/sugrwe.csv') %>% as_tibble() %>% 
  mutate(., vict = number-survivor)

# Full model on whole dataset #

mod.surv <- glm(
  cbind(survivor, vict) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = quasibinomial,
  data = dat.org
)

summary(mod.surv)

mod.surv1 <- glm(
  survival ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  family = quasibinomial,
  data = dat.org
)

summary(mod.surv1)

# Checking means

emm.temp <- emmeans(mod.surv1, ~ temp, type = "response")
emm.swp <- emmeans(mod.surv1, ~ lgswp, type = "response")
emm.cu <- emmeans(mod.surv1, ~ sqrcu, type = "response")

emm.surv <- emmeans(
  mod.surv1,
  ~ temp * lgswp,
  type = "response"
)

plot(emm.surv)

emmip(
  mod.surv1,
  lgswp ~ temp,
  type = "response"
)


emm.swp <- emmeans(
  mod.surv1,
  ~ lgswp | temp,
  type = "response"
)

emm.cu <- emmeans(
  mod.surv1,
  ~ sqrcu | temp,
  type = "response"
)

emmip(mod.surv1,
      lgswp ~ temp,
      type = "response")

emmip(mod.surv1,
      sqrcu ~ temp,
      type = "response")

# comparing simple models

mod.surv0 <- glm(
  survival ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu),
  family = quasibinomial,
  data = dat.org
)

mod.surv1 <- glm(
  survival ~ (factor(temp) +
                factor(lgswp) +
                factor(sqrcu))^2,
  family = quasibinomial,
  data = dat.org
)

anova(mod.surv0, mod.surv1, test = "F")

drop1(mod.surv1, test = "F")

emmip(
  mod.surv1,
  lgswp ~ temp,
  type = "response",
  CIs = TRUE
) +
  theme_bw() +
  labs(
    x = "Temperature (°C)",
    y = "Predicted survival (%)",
    colour = "SWP (kPa)"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

surv_swp <- as.data.frame(
  emmeans(mod.surv1,
          ~ temp * lgswp,
          type = "response")
)

ggplot(
  surv_swp,
  aes(
    x = temp,
    y = prob,
    colour = factor(lgswp),
    group = lgswp
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.1
  ) +
  theme_bw() +
  labs(
    x = "Temperature (°C)",
    y = "Predicted survival",
    colour = "SWP (kPa)"
  )

ggplot(
  surv_swp,
  aes(
    x = temp,
    y = prob,
    group = 1
  )
) +
  geom_line() +
  geom_point(size = 2.5) +
  geom_errorbar(
    aes(
      ymin = asymp.LCL,
      ymax = asymp.UCL
    ),
    width = 0.1
  ) +
  facet_wrap(~ lgswp) +
  theme_bw() +
  labs(
    x = "Temperature (°C)",
    y = "Predicted survival"
  )

###################### separate in cold and warm ##############################

dat.cold <- dat.org %>%
  filter(temp %in% c(20.1, 9.7, 4.5))


dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1, 27.7))

# cold

mod.surv.cold.binom <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = binomial,
  data = dat.cold
)

summary(mod.surv.cold.binom)

# calculate dispersion

dispersion <- sum(residuals(mod.surv.cold.binom, type = "pearson")^2) /
  mod.surv.cold.binom$df.residual

dispersion

# dispersion = 1.67 - mild to moderate dispersion so quasibinomial is better

mod.surv.cold <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = quasibinomial,
  data = dat.cold
)

summary(mod.surv.cold)

drop1(mod.surv.cold, test = "F")


# dispersion = 1.67 indicating moderate overdispersion
# not severe, but better to use quasibinomial error distribution

mod.surv.cold2 <- glm(
  cbind(survivor, number - survivor) ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  family = quasibinomial,
  data = dat.cold
)


mod.surv.cold.final <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) + factor(lgswp) + factor(sqrcu),
  family = quasibinomial,
  data = dat.cold
)

summary(mod.surv.cold.final)
drop1(mod.surv.cold.final, test = "F")

mod.surv.cold2 <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp) +
    factor(temp):factor(sqrcu),
  family = quasibinomial,
  data = dat.cold
)

summary(mod.surv.cold2)


drop1(mod.surv.cold2, test = "F")
drop1(mod.surv.cold1, test = "F")

par(mfrow = c(2, 2))
plot(mod.surv.cold.final)

# check diagnostics plot

library(DHARMa)

sim <- simulateResiduals(mod.surv.cold.final)

plot(sim)

testDispersion(sim)
testOutliers(sim)
testUniformity(sim)

# Overdispersion
deviance(mod.surv.cold.final) / df.residual(mod.surv.cold.final)

par(mfrow = c(1, 2))

# Influence
plot(cooks.distance(mod.surv.cold.final))

# Pearson residuals
plot(
  fitted(mod.surv.cold.final),
  residuals(mod.surv.cold.final, type = "pearson")
)
abline(h = 0, lty = 2)

anova(mod.surv.cold, mod.surv.cold1, mod.surv.cold.final, mod.surv.cold2, test = "F")


########################### Warm ###############################################

mod.surv.warm <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = quasibinomial,
  data = dat.warm
)

summary(mod.surv.warm)

drop1(mod.surv.warm, test = "F")

mod.surv.warm.binom <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = binomial,
  data = dat.warm
)

dispersion <- sum(residuals(mod.surv.warm.binom, type = "pearson")^2) /
  mod.surv.warm.binom$df.residual

dispersion

# dispersion = 2.18 indicating moderate overdispersion
# best to use quasibinomial error distribution

mod.surv.warm1 <- glm(
  cbind(survivor, number - survivor) ~ (factor(temp) + 
                                          factor(lgswp) + 
                                          factor(sqrcu))^2,
  family = quasibinomial,
  data = dat.warm
)

summary(mod.surv.warm1)

mod.surv.warm2 <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) + factor(lgswp) + factor(sqrcu),
  family = quasibinomial,
  data = dat.warm
)

summary(mod.surv.warm2)

anova(mod.surv.warm, mod.surv.warm1, mod.surv.warm2, test = "F")

drop1(mod.surv.warm1, test = "F")

mod.surv.warm.final <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp) +
    factor(temp):factor(sqrcu),
  family = quasibinomial,
  data = dat.warm
)

summary(mod.surv.warm.final)

drop1(mod.surv.warm.final, test = "F")

# Overdispersion
deviance(mod.surv.warm.final) / df.residual(mod.surv.warm.final)

par(mfrow = c(1, 2))

# Influence
plot(cooks.distance(mod.surv.warm.final))

# Pearson residuals
plot(
  fitted(mod.surv.warm.final),
  residuals(mod.surv.warm.final, type = "pearson")
)
abline(h = 0, lty = 2)

anova(mod.surv.warm, mod.surv.warm1, mod.surv.warm.final, mod.surv.warm2, test = "F")

# The best model is the mod.surv.warm.final containing interactions with temperature


