
# figure 1 and 2 in main manuscript - internal Cu and Pb concentrations

library(tidyverse) # v. 2.0.0

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

# Figures for internal Cu and Pb concentrations

dat.cold <- dat.org %>%
  filter(tem %in% c(20.1, 9.7, 4.5))

dat.warm <- dat.org %>%
  filter(tem %in% c(20.1, 24.1, 27.7))

# Important: The statistical analysis has been performed on the log-transformed
# scale to address non-normality issues. However, to improve interpretation in
# figures, the original concentration scale is used.

mod.warm.lm <- lm(
  cu ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.warm
)

mod.cold.lm <- lm(
  cu ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.cold
)

library(emmeans)

emm.temp.warm <- emmeans(mod.warm.lm, ~ tem)
emm.swp.warm <- emmeans(mod.warm.lm, ~ lgswp)
emm.cu.warm <- emmeans(mod.warm.lm, ~ sqrcu)

temp_warm.df <- as.data.frame(emm.temp.warm)
swp_warm.df <- as.data.frame(emm.swp.warm)
cu_warm.df <- as.data.frame(emm.cu.warm)

emm.temp.cold <- emmeans(mod.cold.lm, ~ tem)
emm.swp.cold <- emmeans(mod.cold.lm, ~ lgswp)
emm.cu.cold <- emmeans(mod.cold.lm, ~ sqrcu)

temp_cold.df <- as.data.frame(emm.temp.cold)
swp_cold.df <- as.data.frame(emm.swp.cold)
cu_cold.df <- as.data.frame(emm.cu.cold)

temp_cold.df$Gradient <- "Cold"
temp_warm.df$Gradient <- "Warm"

temp.df <- rbind(temp_cold.df, temp_warm.df)

swp_cold.df$Gradient <- "Cold"
swp_warm.df$Gradient <- "Warm"

swp.df <- rbind(swp_cold.df, swp_warm.df)

cu_cold.df$Gradient <- "Cold"
cu_warm.df$Gradient <- "Warm"

cu.df <- rbind(cu_cold.df, cu_warm.df)


ymax <- max(
  temp.df$upper.CL,
  swp.df$upper.CL,
  cu.df$upper.CL,
  na.rm = TRUE
)

ymin <- min(
  temp.df$lower.CL,
  swp.df$lower.CL,
  cu.df$lower.CL,
  na.rm = TRUE
)

p.temp <- ggplot(temp.df,
                 aes(x = factor(tem),
                     y = emmean,
                     colour = Gradient,
                     group = Gradient)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.1) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  scale_colour_manual(
    values = c("Cold" = "steelblue4",
               "Warm" = "red4")
  ) +
  labs(x = "Temperature (°C)",
       y = "Predicted internal\nCu concentration (mg/kg)") +
  theme_bw()


p.swp <- ggplot(swp.df,
                aes(x = factor(lgswp),
                    y = emmean,
                    colour = Gradient,
                    group = Gradient)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.1) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  scale_colour_manual(
    values = c("Cold" = "steelblue4",
               "Warm" = "red4")
  ) +
  labs(x = "Soil water potential (-kPa)",
       y = "Predicted internal Cu concentration (mg/kg)") +
  theme_bw()

p.swp <- p.swp +
  labs(y = NULL) +
  theme(
    axis.title.y = element_blank()
  )


p.cu <- ggplot(cu.df,
               aes(x = factor(sqrcu),
                   y = emmean,
                   colour = Gradient,
                   group = Gradient)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.1) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  scale_colour_manual(
    values = c("Cold" = "steelblue4",
               "Warm" = "red4")
  ) +
  labs(x = "Soil Cu concentration (mg/kg)",
       y = "Predicted internal Cu concentration") +
  theme_bw()

p.cu <- p.cu +
  labs(y = NULL) +
  theme(
    axis.title.y = element_blank()
  )

common_theme <- theme_bw() +
  theme(
    legend.position = "bottom"
  )

p.temp <- p.temp + common_theme
p.swp <- p.swp + common_theme
p.cu <- p.cu + common_theme

library(patchwork)

p.all <- (p.temp + p.swp + p.cu) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "cu_internal.pdf",
  plot = p.all,
  width = 240,
  height = 105,
  units = "mm"
)
################################################################################
##################### effect size forest plot###################################
################################################################################


################################################################################
############################# cold #############################################

library(broom)

mod.cold.lm <- lm(
  log1p(cu) ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.cold
)

emm.tem <- emmeans(mod.cold.lm, ~ tem)
emm.swp <- emmeans(mod.cold.lm, ~ lgswp)
emm.cu <- emmeans(mod.cold.lm, ~ sqrcu)

# forest plot for supplementary comparing treatments
tem.con <- contrast(
  emm.tem,
  method = "trt.vs.ctrl",
  ref = 3
)

swp.con <- contrast(
  emm.swp,
  method = "trt.vs.ctrl",
  ref = 1 # modify to your SWP reference
)

cu.con <- contrast(
  emm.cu,
  method = "trt.vs.ctrl",
  ref = 1 # first Cu level
)

tem.df <- as.data.frame(summary(tem.con, infer = TRUE))
swp.df <- as.data.frame(summary(swp.con, infer = TRUE))
cu.df <- as.data.frame(summary(cu.con, infer = TRUE))

tem.df$group <- "Temperature"
swp.df$group <- "SWP"
cu.df$group <- "Cu contamination"

coef.df <- rbind(tem.df, swp.df, cu.df)

coef.df$contrast <- recode(
  coef.df$contrast,
  
  "tem4.5 - tem20.1" = "4.5 °C vs 20.1 °C",
  "tem9.7 - tem20.1" = "9.7 °C vs 20.1 °C",
  
  "lgswp5 - lgswp3.2" = "SWP 5 vs 3.2",
  "lgswp10 - lgswp3.2" = "SWP 10 vs 3.2",
  "lgswp39.8 - lgswp3.2" = "SWP 39.8 vs 3.2",
  "lgswp63.1 - lgswp3.2" = "SWP 63.1 vs 3.2",
  "lgswp100 - lgswp3.2" = "SWP 100 vs 3.2",
  
  "sqrcu1629 - sqrcu17" = "High Cu vs Low Cu",
  "sqrcu436 - sqrcu17" = "Medium Cu vs Low Cu"
)

coef.df$sig <- coef.df$lower.CL > 0 |
  coef.df$upper.CL < 0

ggplot(
  coef.df,
  aes(y = reorder(contrast, estimate),
      x = estimate,
      colour = sig)
) +
  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey50") +
  geom_point(size = 4.5) +
  geom_errorbarh(
    aes(xmin = lower.CL,
        xmax = upper.CL),
    height = 0.2,
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "grey50",
      "TRUE" = "blue4"
    ),
    guide = "none"
  ) +
  facet_grid(group ~ .,
             scales = "free_y",
             space = "free_y") +
  labs(
    x = "Difference from reference treatment",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )


################################################################################
######### warm #################################################################


mod.warm.lm <- lm(
  log1p(cu) ~ (factor(tem) + factor(lgswp) + factor(sqrcu))^2,
  data = dat.warm
)

emm.tem <- emmeans(mod.warm.lm, ~ tem)
emm.swp <- emmeans(mod.warm.lm, ~ lgswp)
emm.cu <- emmeans(mod.warm.lm, ~ sqrcu)
#emm.int.temp.swp <- emmeans(
#  mod.warm.lm,
#  ~ lgswp | tem,
#  type = "response"
#)
#emm.int.cu.swp <- emmeans(
#  mod.warm.lm,
#  ~ lgswp | sqrcu,
#  type = "response"
#)

# forest plot for supplementary comparing treatments
tem.con <- contrast(
  emm.tem,
  method = "trt.vs.ctrl",
  ref = 1
)

swp.con <- contrast(
  emm.swp,
  method = "trt.vs.ctrl",
  ref = 1 # modify to your SWP reference
)

cu.con <- contrast(
  emm.cu,
  method = "trt.vs.ctrl",
  ref = 1 # first Cu level
)

#swp.simple <- contrast(
#  emm.int.temp.swp,
#  method = "trt.vs.ctrl",
#  ref = 1
#)

#cu.simple <- contrast(
#  emm.int.cu.swp,
#  method = "trt.vs.ctrl",
#  ref = 1
#)

tem.df <- as.data.frame(summary(tem.con, infer = TRUE))
swp.df <- as.data.frame(summary(swp.con, infer = TRUE))
cu.df <- as.data.frame(summary(cu.con, infer = TRUE))
#temp.swp.simple.df <- as.data.frame(
#  confint(swp.simple)
#)
#cu.swp.simple.df <- as.data.frame(
#  confint(cu.simple)
#)

tem.df$group <- "Temperature"
swp.df$group <- "SWP"
cu.df$group <- "Cu contamination"
#temp.swp.simple.df$group <- paste0(
#  "SWP at ",
#  temp.swp.simple.df$tem,
#  " °C"
#)
#cu.swp.simple.df$group <- paste0(
#  "SWP at ",
#  cu.swp.simple.df$sqrcu,
#  "Cu at"
#)

#coef.df <- rbind(tem.df, swp.df, cu.df, temp.swp.simple.df, cu.swp.simple.df)

coef.df <- rbind(tem.df, swp.df, cu.df)
                 
coef.df$contrast <- recode(
  coef.df$contrast,
  
  "tem27.7 - tem20.1" = "27.7 °C vs 20.1 °C",
  "tem24.1 - tem20.1" = "24.1 °C vs 20.1 °C",
  
  "lgswp5 - lgswp3.2" = "SWP 5 vs 3.2",
  "lgswp10 - lgswp3.2" = "SWP 10 vs 3.2",
  "lgswp39.8 - lgswp3.2" = "SWP 39.8 vs 3.2",
  "lgswp63.1 - lgswp3.2" = "SWP 63.1 vs 3.2",
  "lgswp100 - lgswp3.2" = "SWP 100 vs 3.2",
  
  "sqrcu1629 - sqrcu17" = "High Cu vs Low Cu",
  "sqrcu436 - sqrcu17" = "Medium Cu vs Low Cu"
)


coef.df$sig <- coef.df$lower.CL > 0 |
  coef.df$upper.CL < 0

ggplot(
  coef.df,
  aes(y = reorder(contrast, estimate),
      x = estimate,
      colour = sig)
) +
  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey50") +
  geom_point(size = 4.5) +
  geom_errorbarh(
    aes(xmin = lower.CL,
        xmax = upper.CL),
    height = 0.2,
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "grey50",
      "TRUE" = "red4"
    ),
    guide = "none"
  ) +
  facet_grid(group ~ .,
             scales = "free_y",
             space = "free_y") +
  labs(
    x = "Difference from reference treatment",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )



################################################################################
################## figure for pb ################################################
################################################################################

dat.cold <- dat.org %>%
  filter(tem %in% c(20.1, 9.7, 4.5))

mod.cold.lm <- lm(
 pb ~ factor(tem) + factor(lgswp) + factor(sqrcu),
  data = dat.cold
)

dat.warm <- dat.org %>%
  filter(tem %in% c(20.1, 24.1, 27.7))

mod.warm.lm <- lm(
  pb ~ factor(tem) + factor(lgswp) + factor(sqrcu),
  data = dat.warm
)

# with pb, the model converges.

# now create plots

library(emmeans)

emm.temp.warm <- emmeans(mod.warm.lm, ~ tem)
emm.swp.warm <- emmeans(mod.warm.lm, ~ lgswp)
emm.cu.warm <- emmeans(mod.warm.lm, ~ sqrcu)

temp_warm.df <- as.data.frame(emm.temp.warm)
swp_warm.df <- as.data.frame(emm.swp.warm)
cu_warm.df <- as.data.frame(emm.cu.warm)

emm.temp.cold <- emmeans(mod.cold.lm, ~ tem)
emm.swp.cold <- emmeans(mod.cold.lm, ~ lgswp)
emm.cu.cold <- emmeans(mod.cold.lm, ~ sqrcu)

temp_cold.df <- as.data.frame(emm.temp.cold)
swp_cold.df <- as.data.frame(emm.swp.cold)
cu_cold.df <- as.data.frame(emm.cu.cold)

temp_cold.df$Gradient <- "Cold"
temp_warm.df$Gradient <- "Warm"

temp.df <- rbind(temp_cold.df, temp_warm.df)

swp_cold.df$Gradient <- "Cold"
swp_warm.df$Gradient <- "Warm"

swp.df <- rbind(swp_cold.df, swp_warm.df)

cu_cold.df$Gradient <- "Cold"
cu_warm.df$Gradient <- "Warm"

cu.df <- rbind(cu_cold.df, cu_warm.df)


ymax <- max(
  temp.df$upper.CL,
  swp.df$upper.CL,
  cu.df$upper.CL,
  na.rm = TRUE
)

ymin <- min(
  temp.df$lower.CL,
  swp.df$lower.CL,
  cu.df$lower.CL,
  na.rm = TRUE
)

p.temp <- ggplot(temp.df,
                 aes(x = factor(tem),
                     y = emmean,
                     colour = Gradient,
                     group = Gradient)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.1) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  scale_colour_manual(
    values = c("Cold" = "steelblue4",
               "Warm" = "red4")
  ) +
  labs(x = "Temperature (°C)",
       y = "Predicted internal\nPb concentration (mg/kg)") +
  theme_bw()


p.swp <- ggplot(swp.df,
                aes(x = factor(lgswp),
                    y = emmean,
                    colour = Gradient,
                    group = Gradient)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.1) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  scale_colour_manual(
    values = c("Cold" = "steelblue4",
               "Warm" = "red4")
  ) +
  labs(x = "Soil water potential (-kPa)",
       y = "Predicted internal Cu concentration (mg/kg)") +
  theme_bw()

p.swp <- p.swp +
  labs(y = NULL) +
  theme(
    axis.title.y = element_blank()
  )


p.cu <- ggplot(cu.df,
               aes(x = factor(sqrcu),
                   y = emmean,
                   colour = Gradient,
                   group = Gradient)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = lower.CL,
                    ymax = upper.CL),
                width = 0.1) +
  coord_cartesian(ylim = c(ymin, ymax)) +
  scale_colour_manual(
    values = c("Cold" = "steelblue4",
               "Warm" = "red4")
  ) +
  labs(x = "Soil Cu concentration (mg/kg)",
       y = "Predicted internal Cu concentration") +
  theme_bw()

p.cu <- p.cu +
  labs(y = NULL) +
  theme(
    axis.title.y = element_blank()
  )

common_theme <- theme_bw() +
  theme(
    legend.position = "bottom"
  )

p.temp <- p.temp + common_theme
p.swp <- p.swp + common_theme
p.cu <- p.cu + common_theme

library(patchwork)

p.all <- (p.temp + p.swp + p.cu) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  "pb_internal.pdf",
  plot = p.all,
  width = 240,
  height = 105,
  units = "mm"
)
################################################################################
##################### effect size forest plot###################################
################################################################################

library(broom)

################################################################################
############################# cold #############################################

mod.cold.lm <- lm(
  log1p(pb) ~ factor(tem) + factor(lgswp) + factor(sqrcu),
  data = dat.cold
)

emm.tem <- emmeans(mod.cold.lm, ~ tem)
emm.swp <- emmeans(mod.cold.lm, ~ lgswp)
emm.cu <- emmeans(mod.cold.lm, ~ sqrcu)

# forest plot for supplementary comparing treatments
tem.con <- contrast(
  emm.tem,
  method = "trt.vs.ctrl",
  ref = 3
)

swp.con <- contrast(
  emm.swp,
  method = "trt.vs.ctrl",
  ref = 1 # modify to your SWP reference
)

cu.con <- contrast(
  emm.cu,
  method = "trt.vs.ctrl",
  ref = 1 # first Cu level
)

tem.df <- as.data.frame(summary(tem.con, infer = TRUE))
swp.df <- as.data.frame(summary(swp.con, infer = TRUE))
cu.df <- as.data.frame(summary(cu.con, infer = TRUE))

tem.df$group <- "Temperature"
swp.df$group <- "SWP"
cu.df$group <- "Cu contamination"

coef.df <- rbind(tem.df, swp.df, cu.df)

coef.df$contrast <- dplyr::recode(
  coef.df$contrast,
  
  "tem4.5 - tem20.1" = "4.5 °C vs 20.1 °C",
  "tem9.7 - tem20.1" = "9.7 °C vs 20.1 °C",
  
  "lgswp5 - lgswp3.2" = "SWP 5 vs 3.2",
  "lgswp10 - lgswp3.2" = "SWP 10 vs 3.2",
  "lgswp39.8 - lgswp3.2" = "SWP 39.8 vs 3.2",
  "lgswp63.1 - lgswp3.2" = "SWP 63.1 vs 3.2",
  "lgswp100 - lgswp3.2" = "SWP 100 vs 3.2",
  
  "sqrcu1629 - sqrcu17" = "High Cu vs Low Cu",
  "sqrcu436 - sqrcu17" = "Medium Cu vs Low Cu"
)

coef.df$sig <- coef.df$lower.CL > 0 |
  coef.df$upper.CL < 0

ggplot(
  coef.df,
  aes(y = reorder(contrast, estimate),
      x = estimate,
      colour = sig)
) +
  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey50") +
  geom_point(size = 4.5) +
  geom_errorbarh(
    aes(xmin = lower.CL,
        xmax = upper.CL),
    height = 0.2,
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "grey50",
      "TRUE" = "blue4"
    ),
    guide = "none"
  ) +
  facet_grid(group ~ .,
             scales = "free_y",
             space = "free_y") +
  labs(
    x = "Difference from reference treatment",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )



################################################################################
######### warm #################################################################

mod.warm.lm <- lm(
  log1p(pb) ~ factor(tem) + factor(lgswp) + factor(sqrcu),
  data = dat.warm
)


emm.tem <- emmeans(mod.warm.lm, ~ tem)
emm.swp <- emmeans(mod.warm.lm, ~ lgswp)
emm.cu <- emmeans(mod.warm.lm, ~ sqrcu)

# forest plot for supplementary comparing treatments
tem.con <- contrast(
  emm.tem,
  method = "trt.vs.ctrl",
  ref = 1
)

swp.con <- contrast(
  emm.swp,
  method = "trt.vs.ctrl",
  ref = 1 # modify to your SWP reference
)

cu.con <- contrast(
  emm.cu,
  method = "trt.vs.ctrl",
  ref = 1 # first Cu level
)

tem.df <- as.data.frame(summary(tem.con, infer = TRUE))
swp.df <- as.data.frame(summary(swp.con, infer = TRUE))
cu.df <- as.data.frame(summary(cu.con, infer = TRUE))

tem.df$group <- "Temperature"
swp.df$group <- "SWP"
cu.df$group <- "Cu contamination"

coef.df <- rbind(tem.df, swp.df, cu.df)

coef.df$contrast <- dplyr::recode(
  coef.df$contrast,
  
  "tem27.7 - tem20.1" = "27.7 °C vs 20.1 °C",
  "tem24.1 - tem20.1" = "24.1 °C vs 20.1 °C",
  
  "lgswp5 - lgswp3.2" = "SWP 5 vs 3.2",
  "lgswp10 - lgswp3.2" = "SWP 10 vs 3.2",
  "lgswp39.8 - lgswp3.2" = "SWP 39.8 vs 3.2",
  "lgswp63.1 - lgswp3.2" = "SWP 63.1 vs 3.2",
  "lgswp100 - lgswp3.2" = "SWP 100 vs 3.2",
  
  "sqrcu1629 - sqrcu17" = "High Cu vs Low Cu",
  "sqrcu436 - sqrcu17" = "Medium Cu vs Low Cu"
)

coef.df$sig <- coef.df$lower.CL > 0 |
  coef.df$upper.CL < 0

ggplot(
  coef.df,
  aes(y = reorder(contrast, estimate),
      x = estimate,
      colour = sig)
) +
  geom_vline(xintercept = 0,
             linetype = 2,
             colour = "grey50") +
  geom_point(size = 4.5) +
  geom_errorbarh(
    aes(xmin = lower.CL,
        xmax = upper.CL),
    height = 0.2,
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "grey50",
      "TRUE" = "red4"
    ),
    guide = "none"
  ) +
  facet_grid(group ~ .,
             scales = "free_y",
             space = "free_y") +
  labs(
    x = "Difference from reference treatment",
    y = NULL
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )


