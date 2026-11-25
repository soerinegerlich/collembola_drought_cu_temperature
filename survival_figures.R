################################################################################
# Survival of Folsomia candida: full factorial design
#
# - Response: number of survivors vs dead per replicate, quasibinomial GLM
#   with logit link (overdispersion relative to the binomial).
# - Analyses run separately for the cold (4.5-20.1 °C) and warm
#   (20.1-27.7 °C) gradients.
# - Model selection: from the full factorial model, interaction terms are
#   removed sequentially (least significant first) using single-term deletion
#   F-tests (drop1), until all remaining highest-order terms have P < 0.05.
#   Main effects are always retained.
# - Tests: Type III F-tests (car::Anova) with sum-to-zero contrasts, using the
#   same (deviance-based) dispersion estimate as drop1.
# - EMMs: survival probabilities averaged over the other factors on the logit
#   scale and back-transformed to the probability scale.
################################################################################

library(tidyverse)  # v. 2.0.0
library(car)
library(emmeans)
library(patchwork)
library(DHARMa)

# Sum-to-zero contrasts: must be set BEFORE any model is fitted
options(contrasts = c("contr.sum", "contr.poly"))

# ---- data ---------------------------------------------------------------------

cu.list  <- c(17, 436, 1629)
swp.list <- c(3, 5, 10, 40, 63, 100)

dat.org <- read_csv("data/sugrwe.csv") %>%
  as_tibble() %>%
  mutate(
    lgswp = swp.list[match(lgswp, sort(unique(lgswp)))],  # log10(kPa) -> kPa
    sqrcu = cu.list[match(sqrcu, sort(unique(sqrcu)))],   # sqrt(mg/kg) -> mg/kg
    dead  = number - survivor
  ) %>%
  filter(number > 0)     # replicates with no individuals carry no information

# check the recoding: 6 SWP levels and 3 Cu levels
count(dat.org, lgswp)
count(dat.org, sqrcu)

dat.cold <- dat.org %>%
  filter(temp %in% c(20.1, 9.7, 4.5))

dat.warm <- dat.org %>%
  filter(temp %in% c(20.1, 24.1, 27.7))

# ---- overdispersion check (Figure S5) ----------------------------------------------
# DHARMa does not support quasi-families, so overdispersion is assessed on the
# binomial version of the full model. Dispersion > 1 justifies quasibinomial.

surv.cold.binom <- glm(
  cbind(survivor, dead) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = binomial,
  data = dat.cold
)
surv.warm.binom <- glm(
  cbind(survivor, dead) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = binomial,
  data = dat.warm
)

pearson_dispersion <- function(model) {
  sum(residuals(model, type = "pearson")^2) / df.residual(model)
}

pearson_dispersion(surv.cold.binom)
pearson_dispersion(surv.warm.binom)

sim.cold <- simulateResiduals(surv.cold.binom)
plot(sim.cold)
testDispersion(sim.cold)

sim.warm <- simulateResiduals(surv.warm.binom)
plot(sim.warm)
testDispersion(sim.warm)

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
# Type III table for the final model (Table S7). Partial eta^2 is not defined
# for GLMs, so F, df and P are reported.
report_glm <- function(full, final, label) {
  cat("\n\n====================", label, "====================\n")
  
  cat("\n-- Three-way interaction in the full model (Table S4) --\n")
  print(car::Anova(full, type = 3, test.statistic = "F",
                   error.estimate = "deviance"))
  
  cat("\n-- Final model structure (Table S3) --\n")
  print(formula(final))
  
  cat("\n-- Type III F-tests, final model (Table S7) --\n")
  print(car::Anova(final, type = 3, test.statistic = "F",
                   error.estimate = "deviance"))
  
  cat("\n-- Dispersion parameter, final model --\n")
  print(summary(final)$dispersion)
}

################################################################################
############################ models ############################################
################################################################################

# ---- cold gradient ----

surv.cold.full <- glm(
  cbind(survivor, dead) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = quasibinomial,
  data = dat.cold
)
surv.cold.final <- simplify_interactions(surv.cold.full, test = "F")
report_glm(surv.cold.full, surv.cold.final, "Survival - cold gradient")

# diagnostics: influence and Pearson residuals
par(mfrow = c(1, 2))
plot(cooks.distance(surv.cold.final))
plot(fitted(surv.cold.final), residuals(surv.cold.final, type = "pearson"))
abline(h = 0, lty = 2)

# ---- warm gradient ----

surv.warm.full <- glm(
  cbind(survivor, dead) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = quasibinomial,
  data = dat.warm
)
surv.warm.final <- simplify_interactions(surv.warm.full, test = "F")
report_glm(surv.warm.full, surv.warm.final, "Survival - warm gradient")

# diagnostics: influence and Pearson residuals
par(mfrow = c(1, 2))
plot(cooks.distance(surv.warm.final))
plot(fitted(surv.warm.final), residuals(surv.warm.final, type = "pearson"))
abline(h = 0, lty = 2)

################################################################################
############ back-transformed EMMs (survival probabilities) ####################
################################################################################

# The logit link has no added constant, so type = "response" back-transforms
# correctly. CI column names differ between families (lower.CL / asymp.LCL),
# so they are standardised here.
get_emm_surv <- function(model, term, gradient) {
  emmeans(model, as.formula(paste("~", term)), type = "response") %>%
    as.data.frame() %>%
    rename(any_of(c(emmean   = "prob",
                    lower.CL = "asymp.LCL",
                    upper.CL = "asymp.UCL"))) %>%
    mutate(Gradient = gradient)
}

collect_emms_surv <- function(mod.cold, mod.warm) {
  terms <- c(temp = "temp", lgswp = "lgswp", sqrcu = "sqrcu")
  lapply(terms, function(tm)
    rbind(get_emm_surv(mod.cold, tm, "Cold"),
          get_emm_surv(mod.warm, tm, "Warm")))
}

emms.surv <- collect_emms_surv(surv.cold.final, surv.warm.final)

# Values for the Results text (survival probability, 95% CI)
lapply(emms.surv, function(d)
  d %>%
    dplyr::select(Gradient, 1, emmean, lower.CL, upper.CL) %>%
    mutate(across(c(emmean, lower.CL, upper.CL), ~ round(.x, 3))))

################################################################################
############################ Figure 3 ##########################################
################################################################################

plot_emm_surv <- function(df, xvar, xlab, ylab = NULL) {
  ggplot(df,
         aes(x = factor(.data[[xvar]]),
             y = emmean,
             colour = Gradient,
             group = Gradient)) +
    geom_point(size = 3) +
    geom_line() +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_colour_manual(values = c("Cold" = "steelblue4",
                                   "Warm" = "red4")) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(legend.position = "bottom")
}

p.temp <- plot_emm_surv(emms.surv$temp, "temp",
                        xlab = "Temperature (°C)",
                        ylab = "Estimated survival probability")
p.swp  <- plot_emm_surv(emms.surv$lgswp, "lgswp",
                        xlab = "Soil water potential (-kPa)")
p.cu   <- plot_emm_surv(emms.surv$sqrcu, "sqrcu",
                        xlab = "Soil Cu concentration (mg/kg)")

fig.surv <- (p.temp + p.swp + p.cu) +
  plot_annotation(tag_levels = "A") +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

fig.surv

ggsave("survival.pdf", plot = fig.surv,
       width = 240, height = 105, units = "mm")

################################################################################
############ SI forest plots: contrasts on the logit scale #####################
################################################################################
# Contrasts of EMMs averaged over the other factors (main-effect contrasts,
# not simple effects at reference levels). Reference levels: 20.1 °C,
# -3 kPa, lowest soil Cu.

label_contrast <- function(x) {
  cu.names <- c("17" = "Low Cu", "436" = "Medium Cu", "1629" = "High Cu")
  x %>%
    str_replace("^temp(\\S+) - temp(\\S+)$", "\\1 °C vs \\2 °C") %>%
    str_replace("^lgswp(\\S+) - lgswp(\\S+)$", "-\\1 vs -\\2 kPa") %>%
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
    rename(any_of(c(lower.CL = "asymp.LCL", upper.CL = "asymp.UCL"))) %>%
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
    geom_errorbar(aes(xmin = lower.CL, xmax = upper.CL),
                  width = 0.2, linewidth = 0.8, orientation = "y") +
    scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = colour),
                        guide = "none") +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Difference from reference treatment\n(logit scale)",
         y = NULL) +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14),
          strip.text = element_text(size = 13, face = "bold"),
          panel.grid.minor = element_blank())
}

fp.surv.cold <- forest_plot(surv.cold.final, "blue4")
fp.surv.warm <- forest_plot(surv.warm.final, "red4")

fp.surv.cold
fp.surv.warm
