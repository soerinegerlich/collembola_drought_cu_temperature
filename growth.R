################################################################################
# Growth (body weight) of Folsomia candida: full factorial design
#
# - Response: log-transformed body weight, linear models (lm). A Gamma GLM
#   with log link gave similar diagnostics; lm on log(weight) was chosen as
#   the simpler model.
# - Analyses run separately for the cold (4.5-20.1 °C) and warm
#   (20.1-27.7 °C) gradients.
# - Model selection: from the full factorial model, interaction terms are
#   removed sequentially (least significant first) using single-term deletion
#   F-tests (drop1), until all remaining highest-order terms have P < 0.05.
#   Main effects are always retained.
# - Tests: Type III F-tests (car::Anova) with sum-to-zero contrasts.
# - EMMs: averaged over the other factors on the log scale and
#   back-transformed to mg (geometric means); contrasts are reported as
#   ratios (fold changes).
################################################################################

library(tidyverse)  # v. 2.0.0
library(car)
library(effectsize)
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

# ---- distribution of the raw data -------------------------------------------------

summary(dat.org$growth)

par(mfrow = c(1, 2))
hist(dat.org$growth, breaks = 30, main = "Body weight", xlab = "Body weight")
hist(log(dat.org$growth), breaks = 30, main = "log(body weight)",
     xlab = "log(body weight)")

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
