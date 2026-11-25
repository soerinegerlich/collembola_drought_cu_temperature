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
# - Tests: Type III F-tests (car::Anova) with sum-to-zero contrasts.
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

pearson_dispersion(surv.cold.binom)   # dispersion = 1.67
pearson_dispersion(surv.warm.binom)   # dispersion = 2.18

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
# Type III table for the final model. Partial eta^2 is not defined for GLMs,
# so F, df and P are reported.
report_glm <- function(full, final, label) {
  cat("\n\n====================", label, "====================\n")
  
  cat("\n-- Three-way interaction in the full model (Table S4) --\n")
  print(car::Anova(full, type = 3, test.statistic = "F", error.estimate = "deviance"))
  
  cat("\n-- Final model structure (Table S3) --\n")
  print(formula(final))
  
  cat("\n-- Type III F-tests, final model --\n")
  print(car::Anova(final, type = 3, test.statistic = "F", error.estimate = "deviance"))
  
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

# diagnostics plot

par(mfrow = c(1, 2))

# Influence
plot(cooks.distance(surv.cold.final))

# Pearson residuals
plot(
  fitted(surv.cold.final),
  residuals(surv.cold.final, type = "pearson")
)
abline(h = 0, lty = 2)


# ---- warm gradient ----

surv.warm.full <- glm(
  cbind(survivor, dead) ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = quasibinomial,
  data = dat.warm
)
surv.warm.final <- simplify_interactions(surv.warm.full, test = "F")
report_glm(surv.warm.full, surv.warm.final, "Survival - warm gradient")

# diagnostics plot

par(mfrow = c(1, 2))

# Influence
plot(cooks.distance(surv.warm.final))

# Pearson residuals
plot(
  fitted(surv.warm.final),
  residuals(surv.warm.final, type = "pearson")
)
abline(h = 0, lty = 2)


################################################################################
############ back-transformed EMMs (survival probabilities) ####################
################################################################################

# CI column names differ between families (lower.CL / asymp.LCL), so they are
# standardised here.
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

