# reproduction figures

# Load library
library(tidyverse) # v. 2.0.0
library(patchwork)
library(emmeans)
library(MASS)

# Input origin data
dat.org <- read_csv('data/sugrwe.csv') %>% as_tibble() %>% mutate(., vict = number-survivor)


dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1))

mod.offspring.final <- glm.nb(
  repro ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp) +
    factor(temp):factor(sqrcu),
  data = dat.warm
)

# EMM figure

# Temperature
emm.temp <- emmeans(mod.offspring.final, ~ temp, type = "response")

temp.df <- as.data.frame(emm.temp)

# SWP
emm.swp <- emmeans(mod.offspring.final, ~ lgswp, type = "response")

swp.df <- as.data.frame(emm.swp)
swp.df$swp <- round(10^swp.df$lgswp, 1)

# Cu
emm.cu <- emmeans(mod.offspring.final, ~ sqrcu, type = "response")

cu.df <- as.data.frame(emm.cu)

# back transform cu contamination
cu.df <- cu.df %>%
  mutate(
    cu_soil = c(17, 436, 1629)[match(
      sqrcu,
      c(4.2, 20.9, 40.4)
    )]
  )

cu.df$Gradient <- "Warm"
swp.df$Gradient <- "Warm"
temp.df$Gradient <- "Warm"

ymin <- 0
ymax <- 700

p.temp <- ggplot(
  temp.df,
  aes(x = factor(temp),
      y = response,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3,
             colour = "red4") +
  geom_line() +
  geom_errorbar(
    aes(ymin = asymp.LCL,
        ymax = asymp.UCL),
    width = 0.1,
    colour = "red4"
  ) +
  scale_colour_manual(
    values = c(
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
      y = response,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3,
             colour = "red4") +
  geom_line() +
  geom_errorbar(
    aes(ymin = asymp.LCL,
        ymax = asymp.UCL),
    width = 0.1,
    colour = "red4"
  ) +
  scale_colour_manual(
    values = c(
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
  aes(x = factor(sqrcu),
      y = response,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3,
             colour = "red4") +
  geom_line() +
  geom_errorbar(
    aes(ymin = asymp.LCL,
        ymax = asymp.UCL),
    width = 0.1,
    colour = "red4"
  ) +
  scale_colour_manual(
    values = c(
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
    legend.position = "none"
  )

ggsave(
  "reproduction.pdf",
  plot = p.all,
  width = 255,
  height = 95,
  units = "mm"
)


###############################################################################

# forest plot

library(broom)
library(dplyr)
library(stringr)
library(ggplot2)

# Extract coefficients
coef.df <- tidy(mod.offspring.final, conf.int = TRUE)

# Estimated marginal means
emm.temp <- emmeans(mod.offspring.final, ~ temp)
emm.swp <- emmeans(mod.offspring.final, ~ lgswp)
emm.cu <- emmeans(mod.offspring.final, ~ sqrcu)

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
forest.df$sig <- forest.df$asymp.LCL > 0 |
  forest.df$asymp.UCL < 0

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
      xmin = asymp.LCL,
      xmax = asymp.UCL
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



