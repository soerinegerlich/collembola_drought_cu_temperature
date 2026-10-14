
# Figure for growth in  main manuscript
library(emmeans)
library(tidyverse) # v. 2.0.0
library(patchwork)

# Input origin data
dat.org <- read_csv('data/sugrwe.csv') %>% as_tibble() %>% mutate(., vict = number-survivor)

dat.cold <- dat.org %>%
  filter(temp %in% c(20.1, 9.7, 4.5))


dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1, 27.7))

mod.growth.warm.twoway <- lm(
  growth ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.warm
)

summary(mod.growth.warm.twoway)

mod.growth.cold.final <- lm(
  growth ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp) +
    factor(temp):factor(sqrcu),
  data = dat.cold
)

summary(mod.growth.cold.final)


# Temperature
emm.temp.cold <- emmeans(mod.growth.cold.final, ~ temp, type = "response")
emm.temp.warm <- emmeans(mod.growth.warm.twoway, ~ temp, type = "response")

temp_cold.df <- as.data.frame(emm.temp.cold)
temp_warm.df <- as.data.frame(emm.temp.warm)

# SWP
emm.swp.cold <- emmeans(mod.growth.cold.final, ~ lgswp, type = "response")
emm.swp.warm <- emmeans(mod.growth.warm.twoway, ~ lgswp, type = "response")

swp_cold.df <- as.data.frame(emm.swp.cold)
swp_warm.df <- as.data.frame(emm.swp.warm)

# Cu
emm.cu.cold <- emmeans(mod.growth.cold.final, ~ sqrcu, type = "response")
emm.cu.warm <- emmeans(mod.growth.warm.twoway, ~ sqrcu, type = "response")

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

ymin <- 0
ymax <- 0.12

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
    y = "Predicted body weight (mg)"
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
  "growth.pdf",
  plot = p.all,
  width = 240,
  height = 105,
  units = "mm"
)

##########

library(broom)
library(dplyr)
library(stringr)
library(ggplot2)

mod.growth.warm.twoway <- lm(
  log(growth) ~ (factor(temp) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.warm
)

mod.growth.cold.final <- lm(
  log(growth) ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp) +
    factor(temp):factor(sqrcu),
  data = dat.cold
)


# Extract coefficients
coef.df <- tidy(mod.growth.warm.twoway, conf.int = TRUE)

# Estimated marginal means
emm.temp <- emmeans(mod.growth.warm.twoway, ~ temp)
emm.swp <- emmeans(mod.growth.warm.twoway, ~ lgswp)
emm.cu <- emmeans(mod.growth.warm.twoway, ~ sqrcu)

# Contrasts relative to reference levels
temp.con <- as.data.frame(
  confint(
    contrast(
      emm.temp,
      method = "trt.vs.ctrl",
      ref = 1
    )
  )
)

swp.con <- as.data.frame(
  confint(
    contrast(
      emm.swp,
      method = "trt.vs.ctrl",
      ref = 1
    )
  )
)

cu.con <- as.data.frame(
  confint(
    contrast(
      emm.cu,
      method = "trt.vs.ctrl",
      ref = 1
    )
  )
)

# Add grouping variable
temp.con$group <- "Temperature"
swp.con$group <- "SWP"
cu.con$group <- "Cu contamination"

# Combine dataframes
forest.df <- bind_rows(
  temp.con,
  swp.con,
  cu.con
)

# Significance indicator
forest.df$sig <- forest.df$lower.CL > 0 |
  forest.df$upper.CL < 0

# Optional: relabel contrasts
forest.df$contrast <- dplyr::recode(
  forest.df$contrast,
  
  "temp24.1 - temp20.1" = "24.1°C vs 20.1°C",
  "temp27.7 - temp20.1" = "27.7°C vs 20.1°C",
  
  "temp20.9 - temp4.2" = "436 vs 17 mg Cu kg⁻¹",
  "temp40.4 - temp4.2" = "1629 vs 17 mg Cu kg⁻¹",
  
  # SWP
  "lgswp0.7 - lgswp0.5" = "-5 vs -3.2 kPa",
  "lgswp1 - lgswp0.5" = "-10 vs -3.2 kPa",
  "lgswp1.6 - lgswp0.5" = "-39.8 vs -3.2 kPa",
  "lgswp1.8 - lgswp0.5" = "-63.1 vs -3.2 kPa",
  "lgswp2 - lgswp0.5" = "-100 vs -3.2 kPa",
  
  # Cu
  "sqrcu20.9 - sqrcu4.2" = "Medium vs Low Cu",
  "sqrcu40.4 - sqrcu4.2" = "High vs Low Cu"
)


# Forest plot
ggplot(
  forest.df,
  aes(
    x = estimate,
    y = reorder(contrast, estimate),
    colour = sig
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey50"
  ) +
  geom_point(size = 4) +
  geom_errorbarh(
    aes(
      xmin = lower.CL,
      xmax = upper.CL
    ),
    height = 0.2
  ) +
  facet_grid(
    group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "grey50",
      "TRUE" = "red4"
    ),
    guide = "none"
  ) +
  labs(
    x = "Difference from reference treatment",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold")
  )

# cold

# Extract coefficients
coef.df <- tidy(mod.growth.cold.final, conf.int = TRUE)

# Estimated marginal means
emm.temp <- emmeans(mod.growth.cold.final, ~ temp)
emm.swp <- emmeans(mod.growth.cold.final, ~ lgswp)
emm.cu <- emmeans(mod.growth.cold.final, ~ sqrcu)

# Contrasts relative to reference levels
temp.con <- as.data.frame(
  confint(
    contrast(
      emm.temp,
      method = "trt.vs.ctrl",
      ref = 3
    )
  )
)

swp.con <- as.data.frame(
  confint(
    contrast(
      emm.swp,
      method = "trt.vs.ctrl",
      ref = 1
    )
  )
)

cu.con <- as.data.frame(
  confint(
    contrast(
      emm.cu,
      method = "trt.vs.ctrl",
      ref = 1
    )
  )
)

# Add grouping variable
temp.con$group <- "Temperature"
swp.con$group <- "SWP"
cu.con$group <- "Cu contamination"

# Combine dataframes
forest.df <- bind_rows(
  temp.con,
  swp.con,
  cu.con
)

# Significance indicator
forest.df$sig <- forest.df$lower.CL > 0 |
  forest.df$upper.CL < 0

# Optional: relabel contrasts
forest.df$contrast <- dplyr::recode(
  forest.df$contrast,
  
  "temp4.5 - temp20.1" = "4.5°C vs 20.1°C",
  "temp9.7 - temp20.1" = "9.7°C vs 20.1°C",
  
  # SWP
  "lgswp0.7 - lgswp0.5" = "-5 vs -3.2 kPa",
  "lgswp1 - lgswp0.5" = "-10 vs -3.2 kPa",
  "lgswp1.6 - lgswp0.5" = "-39.8 vs -3.2 kPa",
  "lgswp1.8 - lgswp0.5" = "-63.1 vs -3.2 kPa",
  "lgswp2 - lgswp0.5" = "-100 vs -3.2 kPa",
  
  # Cu
  "sqrcu20.9 - sqrcu4.2" = "Medium vs Low Cu",
  "sqrcu40.4 - sqrcu4.2" = "High vs Low Cu"
)


# Forest plot
ggplot(
  forest.df,
  aes(
    x = estimate,
    y = reorder(contrast, estimate),
    colour = sig
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2,
    colour = "grey50"
  ) +
  geom_point(size = 4) +
  geom_errorbarh(
    aes(
      xmin = lower.CL,
      xmax = upper.CL
    ),
    height = 0.2
  ) +
  facet_grid(
    group ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "grey50",
      "TRUE" = "blue4"
    ),
    guide = "none"
  ) +
  labs(
    x = "Difference from reference treatment",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold")
  )



