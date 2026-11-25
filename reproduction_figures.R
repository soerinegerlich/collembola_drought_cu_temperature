
# Reproduction figure for main manuscript

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
  )


dat.repro <- dat.org %>%
  filter(temp %in% c(20.1, 24.1))

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

################################################################################
############################ figure ############################################
################################################################################

plot_emm_repro <- function(df, xvar, xlab, ylab = NULL, ylim, line = TRUE) {
  p <- ggplot(df,
              aes(x = factor(.data[[xvar]]),
                  y = emmean,
                  group = 1)) +
    geom_point(size = 3, colour = "red4") +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                  width = 0.1, colour = "red4") +
    coord_cartesian(ylim = ylim) +
    labs(x = xlab, y = ylab) +
    theme_bw()
  
  if (line) p <- p + geom_line(colour = "red4")
  p
}


# common y-range across panels, starting at zero
ylim.repro <- c(0, max(unlist(lapply(emms.repro, `[[`, "upper.CL")),
                       na.rm = TRUE))

p.temp <- plot_emm_repro(emms.repro$temp, "temp",
                         xlab = "Temperature (°C)",
                         ylab = "Estimated number of juveniles",
                         ylim = ylim.repro,
                         line = FALSE)
p.swp  <- plot_emm_repro(emms.repro$lgswp, "lgswp",
                         xlab = "Soil water potential (-kPa)",
                         ylim = ylim.repro)
p.cu   <- plot_emm_repro(emms.repro$sqrcu, "sqrcu",
                         xlab = "Soil Cu concentration (mg/kg)",
                         ylim = ylim.repro)

fig.repro <- (p.temp + p.swp + p.cu) +
  plot_annotation(tag_levels = "A")

fig.repro

ggsave("reproduction.pdf", plot = fig.repro,
       width = 245, height = 95, units = "mm")

################################################################################
############ SI forest plot: contrasts on the log scale (Figure S21) ###########
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
    rename(any_of(c(lower.CL = "asymp.LCL", upper.CL = "asymp.UCL"))) %>%
    mutate(contrast = label_contrast(as.character(contrast)),
           group = group)
}

coef.repro <- bind_rows(
  get_contrasts(repro.final, "temp",  "Temperature"),
  get_contrasts(repro.final, "lgswp", "SWP"),
  get_contrasts(repro.final, "sqrcu", "Cu contamination")
) %>%
  mutate(group = factor(group,
                        levels = c("Temperature", "SWP", "Cu contamination")),
         sig   = lower.CL > 0 | upper.CL < 0)

coef.repro   # check the 24.1 vs 20.1 °C contrast against Table S9

fp.repro <- ggplot(coef.repro,
                   aes(y = reorder(contrast, estimate),
                       x = estimate,
                       colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_point(size = 4.5) +
  geom_errorbarh(aes(xmin = lower.CL, xmax = upper.CL),
                 height = 0.2, linewidth = 0.8) +
  scale_colour_manual(values = c("FALSE" = "grey50", "TRUE" = "red4"),
                      guide = "none") +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Difference from reference treatment\n(log scale)",
       y = NULL) +
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 13, face = "bold"),
        panel.grid.minor = element_blank())

fp.repro




