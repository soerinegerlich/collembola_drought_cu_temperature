
#### Growth figures in main manuscript

library(tidyverse)  # v. 2.0.0
library(emmeans)
library(patchwork)

# Sum-to-zero contrasts: must be set BEFORE any model is fitted
options(contrasts = c("contr.sum", "contr.poly"))

# ---- data ---------------------------------------------------------------------

cu.list  <- c(17, 436, 1629)
swp.list <- c(3, 5, 10, 40, 63, 100)

dat.org <- read_csv("data/sugrwe.csv") %>%
  as_tibble() %>%
  mutate(
    lgswp = swp.list[match(lgswp, sort(unique(lgswp)))],  # log10(kPa) -> kPa
    sqrcu = cu.list[match(sqrcu, sort(unique(sqrcu)))]    # sqrt(mg/kg) -> mg/kg
  ) %>%
  filter(!is.na(growth))   # replicates without surviving animals have no weight

# check the recoding: 6 SWP levels and 3 Cu levels
count(dat.org, lgswp)
count(dat.org, sqrcu)

# log() requires positive values
stopifnot(all(dat.org$growth > 0))

dat.cold <- dat.org %>%
  filter(temp %in% c(20.1, 9.7, 4.5))

dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1, 27.7))

# ---- helper functions ---------------------------------------------------------

# Sequential backward simplification of interaction terms.
# test = "F" for lm and quasibinomial glm; test = "Chisq" for glm.nb
simplify_interactions <- function(model, alpha = 0.05, test = "F") {
  repeat {
    d <- drop1(model, test = test)
    d <- d[grepl(":", rownames(d)), , drop = FALSE]   # interaction terms only
    if (nrow(d) == 0) break
    p <- d[[grep("^Pr", names(d))]]
    if (all(p < alpha, na.rm = TRUE)) break
    worst <- rownames(d)[which.max(p)]
    message("Removing ", worst, " (P = ", signif(max(p), 3), ")")
    model <- update(model, as.formula(paste(". ~ . -", worst)))
  }
  model
}

# Prints the tests for Table S4 (three-way), Table S3 (structure) and the
# Type III table with partial eta^2 for the final model (Table S8).
report_models <- function(full, final, label) {
  cat("\n\n====================", label, "====================\n")
  
  cat("\n-- Three-way interaction in the full model (Table S4) --\n")
  a.full <- car::Anova(full, type = 3)
  print(a.full)
  print(effectsize::eta_squared(a.full, partial = TRUE, ci = NULL))
  
  cat("\n-- Final model structure (Table S3) --\n")
  print(formula(final))
  
  cat("\n-- Type III tests, final model --\n")
  a.final <- car::Anova(final, type = 3)
  print(a.final)
  print(effectsize::eta_squared(a.final, partial = TRUE, ci = NULL))
  
  par(mfrow = c(2, 2))
  plot(final, main = label)   # diagnostics (Figure S6)
}

################################################################################
############################ models ############################################
################################################################################

# ---- cold gradient ----

growth.cold.full <- lm(
  log(growth) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  data = dat.cold
)
growth.cold.final <- simplify_interactions(growth.cold.full)
report_models(growth.cold.full, growth.cold.final, "Growth - cold gradient")

# ---- warm gradient ----

growth.warm.full <- lm(
  log(growth) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  data = dat.warm
)
growth.warm.final <- simplify_interactions(growth.warm.full)
report_models(growth.warm.full, growth.warm.final, "Growth - warm gradient")

################################################################################
############ back-transformed EMMs (geometric means) ###########################
################################################################################

# emmeans recognises the log() in the model formula, so type = "response"
# back-transforms means and CIs with exp().
get_emm_growth <- function(model, term, gradient) {
  emmeans(model, as.formula(paste("~", term)), type = "response") %>%
    as.data.frame() %>%
    rename(emmean = response) %>%
    mutate(Gradient = gradient)
}

collect_emms_growth <- function(mod.cold, mod.warm) {
  terms <- c(temp = "temp", lgswp = "lgswp", sqrcu = "sqrcu")
  lapply(terms, function(tm)
    rbind(get_emm_growth(mod.cold, tm, "Cold"),
          get_emm_growth(mod.warm, tm, "Warm")))
}

emms.growth <- collect_emms_growth(growth.cold.final, growth.warm.final)

# Values for the Results text (back-transformed means and 95% CIs)
lapply(emms.growth, function(d)
  d %>%
    dplyr::select(Gradient, 1, emmean, lower.CL, upper.CL) %>%
    mutate(across(c(emmean, lower.CL, upper.CL), ~ signif(.x, 3))))

# Fold changes between temperatures: ratios of back-transformed means
fold_diff <- function(model, term = "temp") {
  summary(pairs(emmeans(model, as.formula(paste("~", term))), reverse = TRUE),
          infer = TRUE, type = "link") %>%
    as.data.frame() %>%
    transmute(contrast,
              fold    = round(exp(estimate), 2),
              lower   = round(exp(lower.CL), 2),
              upper   = round(exp(upper.CL), 2),
              p.value = signif(p.value, 3))
}

fold_diff(growth.cold.final)
fold_diff(growth.warm.final)

################################################################################
############################ figure ############################################
################################################################################

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

# common y-range across panels, starting at zero
ylim.growth <- c(0, max(unlist(lapply(emms.growth, `[[`, "upper.CL")),
                        na.rm = TRUE))

p.temp <- plot_emm(emms.growth$temp, "temp",
                   xlab = "Temperature (°C)",
                   ylab = "Estimated body dry mass (mg)",
                   ylim = ylim.growth)
p.swp  <- plot_emm(emms.growth$lgswp, "lgswp",
                   xlab = "Soil water potential (-kPa)",
                   ylim = ylim.growth)
p.cu   <- plot_emm(emms.growth$sqrcu, "sqrcu",
                   xlab = "Soil Cu concentration (mg/kg)",
                   ylim = ylim.growth)

fig.growth <- (p.temp + p.swp + p.cu) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig.growth

ggsave("growth.pdf", plot = fig.growth,
       width = 240, height = 105, units = "mm")

################################################################################
############ SI forest plots: contrasts on the log scale #######################
################################################################################
# Contrasts of EMMs averaged over the other factors (main-effect contrasts,
# not simple effects at reference levels). Reference levels: 20.1 °C,
# SWP 3.2 kPa, lowest soil Cu.

label_contrast <- function(x) {
  cu.names <- c("17" = "Low Cu", "436" = "Medium Cu", "1629" = "High Cu")
  x %>%
    str_replace("^temp(\\S+) - temp(\\S+)$", "\\1 °C vs \\2 °C") %>%
    str_replace("^lgswp(\\S+) - lgswp(\\S+)$", "SWP \\1 vs \\2") %>%
    str_replace_all("sqrcu(\\d+)", function(m) cu.names[sub("sqrcu", "", m)]) %>%
    str_replace(" - ", " vs ")
}

get_contrasts <- function(model, term, group) {
  emm <- emmeans(model, as.formula(paste("~", term)))
  lev <- as.character(summary(emm)[[term]])
  ref <- if (term == "temp") which(lev == "20.1") else 1
  
  summary(contrast(emm, method = "trt.vs.ctrl", ref = ref),
          infer = TRUE, type = "link") %>%
    as.data.frame() %>%
    mutate(contrast = label_contrast(as.character(contrast)),
           group = group)
}

forest_plot <- function(model, colour) {
  coef.df <- bind_rows(
    get_contrasts(model, "temp",  "Temperature"),
    get_contrasts(model, "lgswp", "SWP"),
    get_contrasts(model, "sqrcu", "Cu contamination")
  ) %>%
    mutate(group = factor(group,
                          levels = c("Temperature", "SWP", "Cu contamination")),
           sig   = lower.CL > 0 | upper.CL < 0)
  
  ggplot(coef.df,
         aes(y = reorder(contrast, estimate),
             x = estimate,
             colour = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
    geom_point(size = 4.5) +
    geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL),
                   height = 0.2, linewidth = 0.8) +
    scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = colour),
                        guide = "none") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Difference from reference treatment\n(log scale)",
         y = NULL) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 13, face = "bold"),
          panel.grid.minor = element_blank())
}

fp.growth.cold <- forest_plot(growth.cold.final, "blue4")
fp.growth.warm <- forest_plot(growth.warm.final, "red4")

fp.growth.cold
fp.growth.warm

