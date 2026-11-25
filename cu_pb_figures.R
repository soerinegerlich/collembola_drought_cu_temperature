################################################################################
# Figures 1-2 (main manuscript) and SI forest plots:
# internal Cu and Pb concentrations
#
# All statistical analyses were performed on log(x + 1)-transformed internal
# concentrations to address non-normality. Figures 1-2 show estimated marginal
# means and 95% CIs from the final models, back-transformed to the original
# scale (mg/kg) as exp(x) - 1; these correspond to geometric means. The SI
# forest plots show contrasts on the log(x + 1) scale, from the same models.
################################################################################

library(tidyverse) # v. 2.0.0
library(emmeans)
library(patchwork)

# ---- data ---------------------------------------------------------------------

cu.list   <- c(17, 436, 1629)
swp.list  <- c(3, 5, 10, 40, 63, 100)
temp.list <- c(4.5, 9.7, 20.1, 24.1, 27.7)

dat.org <- read_csv("data/cuaspb.csv") %>%
  as_tibble() %>%
  mutate(
    sqrcu = case_when(
      cul == "L1" ~ cu.list[1],
      cul == "L2" ~ cu.list[2],
      cul == "L3" ~ cu.list[3]
    ),
    tem = case_when(
      temp == "5"  ~ temp.list[1],
      temp == "10" ~ temp.list[2],
      temp == "20" ~ temp.list[3],
      temp == "24" ~ temp.list[4],
      temp == "28" ~ temp.list[5]
    ),
    lgswp = case_when(
      swp == "3"   ~ swp.list[1],
      swp == "5"   ~ swp.list[2],
      swp == "10"  ~ swp.list[3],
      swp == "40"  ~ swp.list[4],
      swp == "70"  ~ swp.list[5],
      swp == "100" ~ swp.list[6]
    )
  ) %>%
  dplyr::select(tem, sqrcu, lgswp, weight, cu, pb, repli) %>%
  filter(between(cu, 0, 2000) & between(pb, 0, 2000))

dat.cold <- dat.org %>%
  filter(tem %in% c(20.1, 9.7, 4.5))

dat.warm <- dat.org %>%
  filter(tem %in% c(20.1, 24.1, 27.7))

# ---- final models (Table S3) ------------------------------------------------------
# These four models are used for everything below: figures, Results values
# and SI forest plots.

# Internal Cu, cold gradient: main effects, temperature × Cu and SWP × Cu
mod.cu.cold <- lm(
  log(cu + 1) ~ factor(tem) + factor(lgswp) + factor(sqrcu) +
    factor(tem):factor(sqrcu) + factor(lgswp):factor(sqrcu),
  data = dat.cold
)

# Internal Cu, warm gradient: main effects, temperature × SWP and SWP × Cu
mod.cu.warm <- lm(
  log(cu + 1) ~ factor(tem) + factor(lgswp) + factor(sqrcu) +
    factor(tem):factor(lgswp) + factor(lgswp):factor(sqrcu),
  data = dat.warm
)

# Internal Pb, both gradients: main effects only
mod.pb.cold <- lm(
  log(pb + 1) ~ factor(tem) + factor(lgswp) + factor(sqrcu),
  data = dat.cold
)
mod.pb.warm <- lm(
  log(pb + 1) ~ factor(tem) + factor(lgswp) + factor(sqrcu),
  data = dat.warm
)

# ---- back-transformed estimated marginal means ------------------------------------

# EMMs and 95% CIs are calculated on the log(x + 1) scale and back-transformed
# explicitly as exp(x) - 1. This is done by hand rather than with
# type = "response", because recent versions of emmeans detect log(x + 1)
# automatically, and adding make.tran("genlog", 1) on top subtracted the
# constant twice (all means were 1 mg/kg too low).
# The check stops the script if a model was not fitted on the log(x + 1) scale.
get_emm <- function(model, term, gradient) {
  resp <- deparse(formula(model)[[2]])
  if (!grepl("^log\\(.+ \\+ 1\\)$", resp)) {
    stop("Model response is '", resp, "', expected log(x + 1).")
  }
  summary(emmeans(model, as.formula(paste("~", term))), type = "link") %>%
    as.data.frame() %>%
    mutate(across(c(emmean, lower.CL, upper.CL), ~ exp(.x) - 1),
           Gradient = gradient) %>%
    dplyr::select(-SE)   # the SE is on the log(x + 1) scale
}

collect_emms <- function(mod.cold, mod.warm) {
  terms <- c(tem = "tem", lgswp = "lgswp", sqrcu = "sqrcu")
  lapply(terms, function(tm)
    rbind(get_emm(mod.cold, tm, "Cold"),
          get_emm(mod.warm, tm, "Warm")))
}

emms.cu <- collect_emms(mod.cu.cold, mod.cu.warm)
emms.pb <- collect_emms(mod.pb.cold, mod.pb.warm)

# ---- Figures 1-2 ----------------------------------------------------------------

plot_emm <- function(df, xvar, xlab, ylab = NULL, ylim) {
  ggplot(df,
         aes(x = factor(.data[[xvar]]),
             y = emmean,
             colour = Gradient,
             group = Gradient)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
    coord_cartesian(ylim = ylim) +
    scale_colour_manual(values = c("Cold" = "steelblue4",
                                   "Warm" = "red4")) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(legend.position = "bottom")
}

make_figure <- function(emms, metal) {
  # common y-range across the three panels
  ylim <- range(unlist(lapply(emms, function(d) c(d$lower.CL, d$upper.CL))),
                na.rm = TRUE)
  
  p.temp <- plot_emm(emms$tem, "tem",
                     xlab = "Temperature (°C)",
                     ylab = paste0("Estimated internal\n", metal,
                                   " concentration (mg/kg)"),
                     ylim = ylim)
  p.swp  <- plot_emm(emms$lgswp, "lgswp",
                     xlab = "Soil water potential (-kPa)", ylim = ylim)
  p.cu   <- plot_emm(emms$sqrcu, "sqrcu",
                     xlab = "Soil Cu concentration (mg/kg)", ylim = ylim)
  
  (p.temp + p.swp + p.cu) +
    plot_annotation(tag_levels = "A") +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

fig.cu <- make_figure(emms.cu, "Cu")
fig.pb <- make_figure(emms.pb, "Pb")

fig.cu
fig.pb

ggsave("cu_internal.pdf", plot = fig.cu,
       width = 240, height = 105, units = "mm")
ggsave("pb_internal.pdf", plot = fig.pb,
       width = 240, height = 105, units = "mm")

# ---- values for the Results text --------------------------------------------------

# Back-transformed means and 95% CIs (mg/kg)
report_emms <- function(emms) {
  lapply(emms, function(d)
    d %>%
      dplyr::select(Gradient, 1, emmean, lower.CL, upper.CL) %>%
      mutate(across(c(emmean, lower.CL, upper.CL), ~ round(.x, 1))))
}

report_emms(emms.cu)
report_emms(emms.pb)

# Fold differences between temperatures, two versions from the same model:
# - fold_contrast: back-transformed contrast exp(difference on the log(x + 1)
#   scale), with 95% CI and P-value. Strictly the ratio of (concentration + 1).
# - fold_means:    ratio of the back-transformed means (no CI).
# At high concentrations the two are nearly identical; at low concentrations
# (e.g. Pb at 20.1 °C) the + 1 makes a noticeable difference. Use one of them
# consistently in the Abstract, Results and Discussion.
fold_diff <- function(model, term = "tem") {
  ctr <- summary(pairs(emmeans(model, as.formula(paste("~", term)))),
                 infer = TRUE, type = "link") %>%
    as.data.frame()
  
  means <- summary(emmeans(model, as.formula(paste("~", term))),
                   type = "link") %>%
    as.data.frame()
  means$bt <- exp(means$emmean) - 1
  names(means)[1] <- "level"
  
  # level names in the contrast labels, e.g. "tem4.5 - tem20.1"
  lv <- str_match(as.character(ctr$contrast),
                  paste0("^", term, "(\\S+) - ", term, "(\\S+)$"))
  
  ctr %>%
    transmute(contrast,
              fold_contrast = round(exp(estimate), 2),
              lower         = round(exp(lower.CL), 2),
              upper         = round(exp(upper.CL), 2),
              p.value       = signif(p.value, 3),
              fold_means    = round(
                means$bt[match(lv[, 2], as.character(means$level))] /
                  means$bt[match(lv[, 3], as.character(means$level))], 2))
}

fold_diff(mod.cu.cold)
fold_diff(mod.cu.warm)
fold_diff(mod.pb.cold)
fold_diff(mod.pb.warm)

# ---- SI forest plots (log(x + 1) scale) ------------------------------------------

label_contrast <- function(x) {
  cu.names <- c("17" = "Low Cu", "436" = "Medium Cu", "1629" = "High Cu")
  x %>%
    str_replace("^tem(\\S+) - tem(\\S+)$", "\\1 °C vs \\2 °C") %>%
    str_replace("^lgswp(\\S+) - lgswp(\\S+)$", "-\\1 vs -\\2 kPa") %>%
    str_replace_all("sqrcu(\\d+)", function(m) cu.names[sub("sqrcu", "", m)]) %>%
    str_replace(" - ", " vs ")
}

# Treatment-vs-reference contrasts on the log(x + 1) scale, averaged over the
# levels of the other factors. Reference levels: 20.1 °C, -3 kPa, lowest soil Cu.
get_contrasts <- function(model, term, group) {
  emm <- emmeans(model, as.formula(paste("~", term)))
  lev <- as.character(summary(emm)[[term]])
  ref <- if (term == "tem") which(lev == "20.1") else 1
  
  summary(contrast(emm, method = "trt.vs.ctrl", ref = ref),
          infer = TRUE, type = "link") %>%
    as.data.frame() %>%
    mutate(contrast = as.character(contrast),
           group = group)
}

forest_plot <- function(model, colour) {
  coef.df <- bind_rows(
    get_contrasts(model, "tem",   "Temperature"),
    get_contrasts(model, "lgswp", "SWP"),
    get_contrasts(model, "sqrcu", "Cu contamination")
  ) %>%
    mutate(
      contrast = label_contrast(contrast),
      group    = factor(group,
                        levels = c("Temperature", "SWP", "Cu contamination")),
      sig      = lower.CL > 0 | upper.CL < 0
    )
  
  ggplot(coef.df,
         aes(y = reorder(contrast, estimate),
             x = estimate,
             colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 4.5) +
    geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL),
                  width = 0.2, linewidth = 0.8, orientation = "y") +
    scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = colour),
                        guide = "none") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Difference from reference treatment\n(log(x + 1) scale)",
         y = NULL) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 13, face = "bold"),
          panel.grid.minor = element_blank())
}

fp.cu.cold <- forest_plot(mod.cu.cold, "blue4")
fp.cu.warm <- forest_plot(mod.cu.warm, "red4")
fp.pb.cold <- forest_plot(mod.pb.cold, "blue4")
fp.pb.warm <- forest_plot(mod.pb.warm, "red4")

fp.cu.cold
fp.cu.warm
fp.pb.cold
fp.pb.warm
