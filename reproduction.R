# looking at reproduction as the response variable

# Load library
library(tidyverse) # v. 2.0.0
library(patchwork)
library(emmeans)


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

library(MASS)

mod.nb <- glm.nb(
  repro ~ factor(temp) *
    factor(lgswp) *
    factor(sqrcu),
  data = dat.org
)

plot(mod.nb)

# splitting reproduction

# First, ask. Did reproduction occur?

dat.org$reproduced <- as.numeric(dat.org$repro > 0)

glm(
  reproduced ~ factor(temp) *
    factor(lgswp) *
    factor(sqrcu),
  family = binomial,
  data = dat.org
)

# Second, if reproduction occurred, how many offspring were produced?

dat.rep <- subset(dat.org,
                  repro > 0)

glm.nb(
  repro ~ factor(temp) *
    factor(lgswp) *
    factor(sqrcu),
  data = dat.rep
)

# Approach:
# 1. Analyse reproduction success (repro > 0) with a binomial model.
# 2. Analyse offspring number among reproducing individuals with a negative 
# binomial model.
# 3.Present reproduction as occurring only within a limited thermal window 
# centered around 20.1-24.1°C.

dat.org$repro_success <- ifelse(dat.org$repro > 0, 1, 0)

mod.repro.success <- glm(
  repro_success ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  family = binomial,
  data = dat.org
)

drop1(mod.repro.success, test = "Chisq")

mod.repro.success.final <- glm(
  repro_success ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp),
  family = binomial,
  data = dat.org
)

anova(
  mod.repro.success.final,
  mod.repro.success,
  test = "Chisq"
)

drop1(mod.repro.success.final, test = "Chisq")

# reproduction occured within a limited thermal window

emm.temp.rep <- emmeans(
  mod.repro.success.final,
  ~ temp,
  type = "response"
)

summary(emm.temp.rep)

emmeans(
  mod.repro.success.final,
  ~ temp * lgswp,
  type = "response"
)

# analyze the number of offspring among reproducing individuals

mod.offspring <- glm.nb(
  repro ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp),
  data = dat.rep
)

summary(mod.offspring)

drop1(mod.offspring, test = "Chisq")

# Cu contamination has no effect. However, we want to keep it

dat.cold <- dat.org %>%
  filter(temp %in% c(20.1, 9.7, 4.5))


dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1, 27.7))


mod.offspring.final <- glm.nb(
  repro ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp),
  data = dat.rep
)

mod.offspring.final <- glm.nb(
  repro ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp),
  data = dat.org
)

summary(mod.offspring.final)

mod.offspring.final.warm <- lm(
  repro ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp),
  data = dat.warm
)

mod.offspring.final.cold <- lm(
  repro ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp),
  data = dat.cold
)

# EMM figure

# Temperature
emm.temp.cold <- emmeans(mod.offspring.final.cold, ~ temp, type = "response")
emm.temp.warm <- emmeans(mod.offspring.final.warm, ~ temp, type = "response")

temp_cold.df <- as.data.frame(emm.temp.cold)
temp_warm.df <- as.data.frame(emm.temp.warm)

# SWP
emm.swp.cold <- emmeans(mod.offspring.final.cold, ~ lgswp, type = "response")
emm.swp.warm <- emmeans(mod.offspring.final.warm, ~ lgswp, type = "response")

swp_cold.df <- as.data.frame(emm.swp.cold)
swp_warm.df <- as.data.frame(emm.swp.warm)

# Cu
emm.cu.cold <- emmeans(mod.offspring.final.cold, ~ sqrcu, type = "response")
emm.cu.warm <- emmeans(mod.offspring.final.warm, ~ sqrcu, type = "response")

cu_cold.df <- as.data.frame(emm.cu.cold)
cu_warm.df <- as.data.frame(emm.cu.warm)

temp_cold.df$Gradient <- "Cold"
temp_warm.df$Gradient <- "Warm"

temp.df <- rbind(temp_cold.df, temp_warm.df)

swp_cold.df$Gradient <- "Cold"
swp_warm.df$Gradient <- "Warm"

swp.df <- rbind(swp_cold.df, swp_warm.df)
swp.df$swp <- round(10^swp.df$lgswp, 1)

cu_cold.df$Gradient <- "Cold"
cu_warm.df$Gradient <- "Warm"

cu.df <- rbind(cu_cold.df, cu_warm.df)
# back transform cu contamination
cu.df <- cu.df %>%
  mutate(
    cu_soil = c(17, 436, 1629)[match(
      sqrcu,
      c(4.2, 20.9, 40.4)
    )]
  )

ymin <- -10
ymax <- 400

p.temp <- ggplot(
  temp.df,
  aes(x = factor(temp),
      y = emmean,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(
    aes(ymin = lower.CL,
        ymax = upper.CL),
    width = 0.1
  ) +
  scale_colour_manual(
    values = c(
      "Cold" = "steelblue4",
      "Warm" = "red4"
    )
  ) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  labs(
    x = "Temperature (°C)",
    y = "Predicted number of offspring"
  ) +
  theme_bw()


p.swp <- ggplot(
  swp.df,
  aes(x = factor(swp),
      y = emmean,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(
    aes(ymin = lower.CL,
        ymax = upper.CL),
    width = 0.1
  ) +
  scale_colour_manual(
    values = c(
      "Cold" = "steelblue4",
      "Warm" = "red4"
    )
  ) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  labs(
    x = "Soil water potential (-kPa)",
    y = NULL
  ) +
  theme_bw()


p.cu <- ggplot(
  cu.df,
  aes(x = factor(cu_soil),
      y = emmean,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(
    aes(ymin = lower.CL,
        ymax = upper.CL),
    width = 0.1
  ) +
  scale_colour_manual(
    values = c(
      "Cold" = "steelblue4",
      "Warm" = "red4"
    )
  ) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  labs(
    x = "Soil Cu concentration (mg/kg)",
    y = NULL
  ) +
  theme_bw()

library(patchwork)

p.all <- (p.temp + p.swp + p.cu) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom"
  )

ggsave(
  "emmeans_reproduction.pdf",
  plot = p.all,
  width = 240,
  height = 105,
  units = "mm"
)

