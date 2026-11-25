################################################################################
# Reproduction of Folsomia candida: full factorial design
#
# - Response: number of juveniles per replicate, negative binomial GLM with
#   log link (strong overdispersion relative to the Poisson).
# - Only 20.1 and 24.1 °C are analysed: egg development takes more than four
#   weeks at 4.5 and 9.7 °C, and reproduction was rare at 27.7 °C.
# - Model selection: from the full factorial model, interaction terms are
#   removed sequentially (least significant first) using single-term deletion
#   likelihood-ratio tests (drop1, test = "Chisq"), until all remaining
#   highest-order terms have P < 0.05. Main effects are always retained.
# - Tests: Type III likelihood-ratio tests (car::Anova) with sum-to-zero
#   contrasts.
# - EMMs: averaged over the other factors on the log scale and
#   back-transformed to numbers of juveniles; contrasts are reported as
#   ratios.
################################################################################

library(tidyverse)  # v. 2.0.0
library(MASS)
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
    sqrcu = cu.list[match(sqrcu, sort(unique(sqrcu)))]    # sqrt(mg/kg) -> mg/kg
  )

# check the recoding: 6 SWP levels and 3 Cu levels
count(dat.org, lgswp)
count(dat.org, sqrcu)

# ---- why only 20.1 and 24.1 °C ----------------------------------------------------

# Mean reproduction and proportion of zeros per temperature
dat.org %>%
  group_by(temp) %>%
  summarise(mean_repro = mean(repro, na.rm = TRUE),
            prop_zero  = mean(repro == 0, na.rm = TRUE))

dat.repro <- dat.org %>%
  filter(temp %in% c(20.1, 24.1))

# ---- choice of error distribution ---------------------------------------------------

# Poisson: dispersion far above 1 -> Poisson inappropriate
repro.pois <- glm(
  repro ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  family = poisson,
  data = dat.repro
)
sum(residuals(repro.pois, type = "pearson")^2) / df.residual(repro.pois)

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

################################################################################
############################ model #############################################
################################################################################

repro.full <- glm.nb(
  repro ~ factor(temp) * factor(lgswp) * factor(sqrcu),
  data = dat.repro
)

repro.final <- simplify_interactions(repro.full, test = "Chisq")

cat("\n-- Three-way interaction in the full model (Table S4) --\n")
car::Anova(repro.full, type = 3, test.statistic = "LR")

cat("\n-- Final model structure (Table S3) --\n")
formula(repro.final)

cat("\n-- Type III likelihood-ratio tests, final model (Table S9) --\n")
car::Anova(repro.final, type = 3, test.statistic = "LR")

cat("\n-- Negative binomial dispersion parameter (theta) --\n")
repro.final$theta

# ---- residual diagnostics (Figure S7) ----------------------------------------------

sim.repro <- simulateResiduals(fittedModel = repro.final, n = 1000)
plot(sim.repro)

testDispersion(sim.repro)
testZeroInflation(sim.repro)
testUniformity(sim.repro)
testOutliers(sim.repro)

################################################################################
############ back-transformed EMMs (numbers of juveniles) ######################
################################################################################

get_emm_repro <- function(model, term) {
  emmeans(model, as.formula(paste("~", term)), type = "response") %>%
    as.data.frame() %>%
    rename(any_of(c(emmean   = "response",
                    lower.CL = "asymp.LCL",
                    upper.CL = "asymp.UCL")))
}

emms.repro <- list(
  temp  = get_emm_repro(repro.final, "temp"),
  lgswp = get_emm_repro(repro.final, "lgswp"),
  sqrcu = get_emm_repro(repro.final, "sqrcu")
)

# Values for the Results text (juveniles per replicate, 95% CI)
lapply(emms.repro, function(d)
  d %>%
    dplyr::select(1, emmean, lower.CL, upper.CL) %>%
    mutate(across(c(emmean, lower.CL, upper.CL), ~ round(.x, 1))))

# Ratio 24.1 vs 20.1 °C (back-transformed contrast)
summary(pairs(emmeans(repro.final, ~ temp), reverse = TRUE),
        infer = TRUE, type = "response")



