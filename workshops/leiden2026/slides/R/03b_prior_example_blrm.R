here::i_am("slides/R/03b_prior_example_blrm.R")

library(here)
library(knitr)
library(tidyverse)
library(ggrepel)
library(latex2exp)
library(patchwork)
library(glue)
library(RBesT)
library(rstan)
library(brms)
library(posterior)
library(tidybayes)
library(bayesplot)
library(gt)
library(ggdist)
library(distributional)
library(mvtnorm)
library(dqrng)
library(emmeans)
library(simsurv)

## ----blrm_setup---------------------------------------------------------------
#| include: false
#| cache: false
# chunks below depend on these definitions, which knitr's cache does not track
blrm_fig_width <- knitr::opts_chunk$get("fig.width")
# half-width columns: 6in figures keep the text readable
knitr::opts_chunk$set(cache = FALSE, fig.width = 6)
dose_grid <- c(0.5, 1, 2, 4, 8, 12, 16)
starting_dose <- 1
ref_dose <- 4
cohort_doses <- starting_dose * 2^(-1:1)
cohort_label <- paste(cohort_doses, collapse = ", ")
cohort_data <- tibble(Dose = cohort_doses, num_patients = 3, num_toxicities = 0)
curve_doses <- seq(0, max(dose_grid), length.out = 101)

# prior on (log alpha, log beta): means m, SDs s, correlation r
blrm_prior <- list(m = c(logit_scaled(0.17), log(1)), s = c(1.7, log(4) / 1.96), r = 0.6)

draw_prior <- function(prior, n = 4000, seed = 4711) {
  S <- diag(prior$s) %*% matrix(c(1, prior$r, prior$r, 1), 2) %*% diag(prior$s)
  x <- withr::with_seed(seed, rmvnorm(n, prior$m, S))
  tibble(.draw = seq_len(n), par1 = x[, 1], log_beta = x[, 2])
}

classify_slope <- function(log_beta) {
  factor(case_when(exp(log_beta) < 1 / 4 ~ "Flat",
                   exp(log_beta) > 4 ~ "Step function",
                   .default = "Reasonable"),
         levels = c("Step function", "Flat", "Reasonable"))
}

# par1 is log(alpha)
blrm_curves <- function(draws, dref = ref_dose, doses = curve_doses) {
  expand_grid(draws, Dose = doses) |>
    mutate(eta = par1 + exp(log_beta) * log(Dose / dref),
           risk = plogis(eta),
           curve_type = classify_slope(log_beta))
}

curve_pct <- function(curves, type) {
  slopes <- distinct(curves, .draw, curve_type)
  round(100 * mean(slopes$curve_type == type), 1)
}

dose_axis <- scale_x_continuous("Dose", breaks = dose_grid,
                                labels = scales::label_number(drop0trailing = TRUE),
                                guide = guide_axis(n.dodge = 2))
risk_axis <- scale_y_continuous("P(DLT)", breaks = c(0, 0.16, 0.33, 0.5, 0.75, 1),
                                labels = scales::percent, limits = c(0, 1))

plot_curve_dist <- function(curves, n_draws = 1000) {
  filter(curves, .draw <= n_draws) |>
    ggplot(aes(Dose, risk, group = .draw, colour = curve_type)) +
    geom_line(alpha = 0.1, linewidth = 0.8) +
    hline_at(c(0.16, 0.33), linetype = 2) +
    scale_colour_manual(values = c("Step function" = "red", "Flat" = "orange",
                                   "Reasonable" = "darkgreen"), drop = FALSE) +
    dose_axis + risk_axis +
    theme_bw(base_size = 18) +
    theme(legend.position = "none")
}

# ribbons, median and EWOC shading (doses with P(overdose) >= 25% are greyed out)
plot_toxicity_curve <- function(curves) {
  s <- curves |>
    summarise(q025 = quantile(risk, 0.025), q25 = quantile(risk, 0.25),
              q50 = median(risk), q75 = quantile(risk, 0.75),
              q975 = quantile(risk, 0.975), p_over = mean(risk > 0.33), .by = Dose)
  ggplot(s, aes(Dose)) +
    geom_rect(data = filter(s, p_over >= 0.25),
              aes(xmin = Dose, xmax = Dose + diff(curve_doses)[1], ymin = 0, ymax = 1),
              fill = "grey60", alpha = 0.5, inherit.aes = FALSE) +
    geom_ribbon(aes(ymin = q025, ymax = q975), fill = "#d1e1ec") +
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "#6497b1") +
    geom_line(aes(y = q50), colour = "#011f4b", linewidth = 1) +
    hline_at(c(0.16, 0.33), linetype = 2) +
    dose_axis + risk_axis +
    theme_bw(base_size = 18)
}

mvn_stanvars <- function(prior) {
  stanvar(prior$m, "M") + stanvar(prior$s, "S") + stanvar(prior$r, "rho") +
    stanvar(scode = "matrix[2, 2] L = [[S[1], 0], [rho * S[2], sqrt(1 - rho^2) * S[2]]];",
            block = "tdata")
}

blrm_model <- bf(num_toxicities | trials(num_patients) ~ logalpha + exp(logbeta) * log(Dose / dref),
                 logalpha + logbeta ~ 1, nl = TRUE, family = binomial())

fit_mvn <- function(model, pars, data, prior, seed = 4711) {
  mvn <- prior_string(paste0("target += multi_normal_cholesky_lpdf([b_", pars[1], "[1], b_",
                             pars[2], "[1]]' | M, L)"), check = FALSE)
  fit <- brm(model, data = mutate(data, dref = ref_dose),
             prior = mvn, stanvars = mvn_stanvars(prior), seed = seed, refresh = 0, silent = 2)
  as_draws_df(fit) |> as_tibble() |>
    transmute(.draw, par1 = .data[[paste0("b_", pars[1], "_Intercept")]],
              log_beta = b_logbeta_Intercept)
}

data_arrows <- list(
  annotate("segment", x = cohort_doses, xend = cohort_doses, y = 0.6, yend = 0.2,
           arrow = arrow(length = unit(0.4, "cm"), type = "closed"), linewidth = 1.5),
  annotate("text", x = mean(range(cohort_doses)) + 1, y = 0.75, size = 6,
           label = paste0("0/3 DLTs\nat ", cohort_label))
)

curves_default <- blrm_curves(draw_prior(blrm_prior))


## ----blrm_prior_curve---------------------------------------------------------
#| echo: false
#| fig-asp: 0.8
plot_toxicity_curve(curves_default)


## ----blrm_prior_intercept-----------------------------------------------------
#| echo: false
#| fig-asp: 0.8
plot_toxicity_curve(curves_default) +
  stat_eye(data = tibble(x = ref_dose, y = plogis(draw_prior(blrm_prior)$par1)),
           aes(x = x, y = y), .width = c(0.5, 0.95), scale = 1, alpha = 0.8,
           inherit.aes = FALSE)


## ----blrm_prior_slope---------------------------------------------------------
#| echo: false
#| fig-asp: 0.8
median_curve <- function(log_slope_mean) {
  blrm_curves(draw_prior(modifyList(blrm_prior, list(m = c(blrm_prior$m[1], log_slope_mean))))) |>
    summarise(risk = median(risk), .by = Dose)
}
plot_toxicity_curve(curves_default) +
  geom_vline(xintercept = ref_dose, linetype = "dashed") +
  geom_line(aes(y = risk), data = median_curve(log(2)), colour = "darkred", linewidth = 2) +
  geom_line(aes(y = risk), data = median_curve(log(1 / 2)), colour = "darkorange", linewidth = 2)


## ----blrm_prior_extreme-------------------------------------------------------
#| echo: false
#| fig-asp: 0.8
plot_toxicity_curve(curves_default) +
  geom_vline(xintercept = ref_dose, linetype = "dashed") +
  geom_line(aes(y = risk), data = median_curve(log(10)), colour = "darkred", linewidth = 2) +
  geom_line(aes(y = risk), data = median_curve(log(1 / 10)), colour = "darkorange", linewidth = 2)


## ----blrm_prior_spaghetti-----------------------------------------------------
#| echo: false
#| fig-asp: 0.8
plot_curve_dist(curves_default)


## ----blrm_prior_wide----------------------------------------------------------
#| echo: false
#| fig-asp: 0.8
wide_prior <- modifyList(blrm_prior, list(s = c(1.7, 3)))
curves_wide <- blrm_curves(draw_prior(wide_prior))
plot_curve_dist(curves_wide)


## ----blrm_post_wide-----------------------------------------------------------
#| echo: false
#| fig-asp: 0.8
curves_wide_post <- blrm_curves(fit_mvn(blrm_model, c("logalpha", "logbeta"), cohort_data, wide_prior))
plot_curve_dist(curves_wide_post) + data_arrows


## ----blrm_post_default--------------------------------------------------------
#| echo: false
#| fig-asp: 0.8
curves_default_post <- blrm_curves(fit_mvn(blrm_model, c("logalpha", "logbeta"), cohort_data, blrm_prior))
plot_curve_dist(curves_default_post) + data_arrows


## ----blrm_brms_code-----------------------------------------------------------
#| eval: false
# blrm_model <- bf(num_toxicities | trials(num_patients) ~ logalpha + exp(logbeta) * log(Dose / dref),
#                  logalpha + logbeta ~ 1, nl = TRUE, family = binomial())
# 
# # brms has no bivariate normal prior across non-linear parameters, so we add it to the target
# prior <- list(m = c(logit_scaled(0.17), log(1)), s = c(1.7, log(4) / 1.96), r = 0.6)
# stanvars <- stanvar(prior$m, "M") + stanvar(prior$s, "S") + stanvar(prior$r, "rho") +
#   stanvar(scode = "matrix[2, 2] L = [[S[1], 0], [rho * S[2], sqrt(1 - rho^2) * S[2]]];",
#           block = "tdata")
# mvn_prior <- prior("target += multi_normal_cholesky_lpdf([b_logalpha[1], b_logbeta[1]]' | M, L)",
#                    check = FALSE)
# 
# data <- tibble(Dose = c(0.5, 1, 2), num_patients = 3, num_toxicities = 0, dref = 4)
# fit <- brm(blrm_model, data = data, prior = mvn_prior, stanvars = stanvars)


## ----blrm_restore_cache-------------------------------------------------------
#| include: false
knitr::opts_chunk$set(cache = TRUE, fig.width = blrm_fig_width)

