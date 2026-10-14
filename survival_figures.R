
# Figure for survival in  main manuscript
library(emmeans)
library(tidyverse) # v. 2.0.0
library(patchwork)

# Input origin data
dat.org <- read_csv('data/sugrwe.csv') %>% as_tibble() %>% 
  mutate(., vict = number-survivor)

dat.cold <- dat.org %>%
  filter(temp %in% c(20.1, 9.7, 4.5))


dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1, 27.7))

# models

# cold
mod.surv.cold2 <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) + factor(lgswp) + factor(sqrcu),
  family = quasibinomial,
  data = dat.cold
)

# warm
mod.surv.warm.final <- glm(
  cbind(survivor, number - survivor) ~ factor(temp) +
    factor(lgswp) +
    factor(sqrcu) +
    factor(temp):factor(lgswp) +
    factor(temp):factor(sqrcu),
  family = quasibinomial,
  data = dat.warm
)

# Temperature
emm.temp.cold <- emmeans(mod.surv.cold2, ~ temp, type = "response")
emm.temp.warm <- emmeans(mod.surv.warm.final, ~ temp, type = "response")

temp_cold.df <- as.data.frame(emm.temp.cold)
temp_warm.df <- as.data.frame(emm.temp.warm)

# SWP
emm.swp.cold <- emmeans(mod.surv.cold2, ~ lgswp, type = "response")
emm.swp.warm <- emmeans(mod.surv.warm.final, ~ lgswp, type = "response")

swp_cold.df <- as.data.frame(emm.swp.cold)
swp_warm.df <- as.data.frame(emm.swp.warm)

# Cu
emm.cu.cold <- emmeans(mod.surv.cold2, ~ sqrcu, type = "response")
emm.cu.warm <- emmeans(mod.surv.warm.final, ~ sqrcu, type = "response")

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

ymin <- 0.5
ymax <- 1.1

p.temp <- ggplot(
  temp.df,
  aes(x = factor(temp),
      y = prob,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(
    aes(ymin = asymp.LCL,
        ymax = asymp.UCL),
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
    y = "Predicted survival probability"
  ) +
  theme_bw()

p.swp <- ggplot(
  swp.df,
  aes(x = factor(swp),
      y = prob,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(
    aes(ymin = asymp.LCL,
        ymax = asymp.UCL),
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
      y = prob,
      colour = Gradient,
      group = Gradient)
) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(
    aes(ymin = asymp.LCL,
        ymax = asymp.UCL),
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
  "surv.pdf",
  plot = p.all,
  width = 240,
  height = 105,
  units = "mm"
)


# create forest plot with pairwise comparison of EMMs
# Contrasts relative to reference levels

# Reference = first row (4.5°C for cold, 20.1°C for warm)
temp_cold.ref <- temp_cold.df$prob[3]
temp_warm.ref <- temp_warm.df$prob[1]

temp_cold.forest <- temp_cold.df[-3, ] %>%
  mutate(
    estimate = prob - temp_cold.ref,
    lower.CL = asymp.LCL - temp_cold.ref,
    upper.CL = asymp.UCL - temp_cold.ref,
    contrast = paste0(temp, " vs ", temp_cold.df$temp[3]),
    group = "Temperature"
  )

temp_warm.forest <- temp_warm.df[-1, ] %>%
  mutate(
    estimate = prob - temp_warm.ref,
    lower.CL = asymp.LCL - temp_warm.ref,
    upper.CL = asymp.UCL - temp_warm.ref,
    contrast = paste0(temp, " vs ", temp_warm.df$temp[1]),
    group = "Temperature"
  )

swp_cold.ref <- swp_cold.df$prob[1]
swp_warm.ref <- swp_warm.df$prob[1]

swp_cold.forest <- swp_cold.df[-1, ] %>%
  mutate(
    estimate = prob - swp_cold.ref,
    lower.CL = asymp.LCL - swp_cold.ref,
    upper.CL = asymp.UCL - swp_cold.ref,
    contrast = paste0(lgswp, " vs ", swp_cold.df$lgswp[1]),
    group = "SWP"
  )

swp_warm.forest <- swp_warm.df[-1, ] %>%
  mutate(
    estimate = prob - swp_warm.ref,
    lower.CL = asymp.LCL - swp_warm.ref,
    upper.CL = asymp.UCL - swp_warm.ref,
    contrast = paste0(lgswp, " vs ", swp_warm.df$lgswp[1]),
    group = "SWP"
  )

cu_cold.ref <- cu_cold.df$prob[1]
cu_warm.ref <- cu_warm.df$prob[1]

cu_cold.forest <- cu_cold.df[-1, ] %>%
  mutate(
    estimate = prob - cu_cold.ref,
    lower.CL = asymp.LCL - cu_cold.ref,
    upper.CL = asymp.UCL - cu_cold.ref,
    contrast = paste0(sqrcu, " vs ", cu_cold.df$sqrcu[1]),
    group = "Cu contamination"
  )

cu_warm.forest <- cu_warm.df[-1, ] %>%
  mutate(
    estimate = prob - cu_warm.ref,
    lower.CL = asymp.LCL - cu_warm.ref,
    upper.CL = asymp.UCL - cu_warm.ref,
    contrast = paste0(sqrcu, " vs ", cu_warm.df$sqrcu[1]),
    group = "Cu contamination"
  )

# Combine ----------------------------------------------------------------

forest.warm.df <- bind_rows(
  temp_warm.forest,
  swp_warm.forest,
  cu_warm.forest
)

forest.warm.df$sig <- forest.warm.df$lower.CL > 0 |
  forest.warm.df$upper.CL < 0

forest.warm.df$contrast <- dplyr::recode(
  forest.warm.df$contrast,
  
  # Temperature
  "24.1 vs 20.1" = "24.1 °C vs 20.1 °C",
  "27.7 vs 20.1" = "27.7 °C vs 20.1 °C",
  
  # SWP
  "0.7 vs 0.5" = "-5 vs -3.2 kPa",
  "1 vs 0.5" = "-10 vs -3.2 kPa",
  "1.6 vs 0.5" = "-39.8 vs -3.2 kPa",
  "1.8 vs 0.5" = "-63.1 vs -3.2 kPa",
  "2 vs 0.5" = "-100 vs -3.2 kPa",
  
  # Cu contamination
  "20.9 vs 4.2" = "Medium vs. Low Cu",
  "40.4 vs 4.2" = "High vs. Low Cu"
)


# Forest plot
ggplot(
  forest.warm.df,
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
  geom_point(size = 3) +
  geom_errorbarh(
    aes(
      xmin = lower.CL,
      xmax = upper.CL
    ),
    height = 0.2,
    linewidth = 0.8
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
    x = "Difference in survival probability compared to reference treatment",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )


############################ cold ##############################################
forest.cold.df <- bind_rows(
  temp_cold.forest,
  swp_cold.forest,
  cu_cold.forest
)

forest.cold.df$sig <- forest.cold.df$lower.CL > 0 |
  forest.cold.df$upper.CL < 0

forest.cold.df$contrast <- dplyr::recode(
  forest.cold.df$contrast,
  
  # Temperature
  "4.5 vs 20.1" = "4.5 °C vs 20.1 °C",
  "9.7 vs 20.1" = "9.7 °C vs 20.1 °C",
  
  # SWP
  "0.7 vs 0.5" = "-5 vs -3.2 kPa",
  "1 vs 0.5" = "-10 vs -3.2 kPa",
  "1.6 vs 0.5" = "-39.8 vs -3.2 kPa",
  "1.8 vs 0.5" = "-63.1 vs -3.2 kPa",
  "2 vs 0.5" = "-100 vs -3.2 kPa",
  
  # Cu contamination
  "20.9 vs 4.2" = "Medium vs. Low Cu",
  "40.4 vs 4.2" = "High vs. Low Cu"
)


# Forest plot
ggplot(
  forest.cold.df,
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
  geom_point(size = 3) +
  geom_errorbarh(
    aes(
      xmin = lower.CL,
      xmax = upper.CL
    ),
    height = 0.2,
    linewidth = 0.8
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
    x = "Difference in survival probability compared to reference treatment",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )

