here::i_am("slides/R/04_dose_finding.R")

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

## -----------------------------------------------------------------------------
#| echo: false
#| warning: false
#| fig.height: 6.5
drfdata <- tibble(
  Dose = seq(0, 500, length.out = 101),
  Emax = Dose / (Dose + 70),
  `sigmoid Emax` = 1 / (1 + (100 / Dose)^3),
  `non-monotone` = 1 /
    (1 + (10 / Dose)^0.9) -
    (Dose > 300) * (Dose - 300)^2 / 120000
) |>
  gather(value = delta, key = Function, -Dose) |>
  mutate(
    pointlabel = ifelse(
      (Dose == 300 & Function == "Emax") |
        (Dose == 100 & Function == "non-monotone") |
        (Dose == 120 & Function == "sigmoid Emax"),
      Function,
      NA
    )
  )

p_base <- drfdata |>
  filter(Function == "Emax") |>
  ggplot(aes(x = Dose, y = delta, col = Function, label = pointlabel)) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9")) +
  theme_bw(base_size = 24) +
  ylab("Percentage of full treatment effect") +
  geom_hline(yintercept = 0, size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(labels = scales::percent) +
  geom_line(size = 3) +
  guides(color = FALSE)

p1 <- p_base +
  geom_text(data = filter(drfdata, Function == "Emax" & !is.na(pointlabel)),
            size = 10, nudge_y=c(-0.1), nudge_x=c(10))

p2 <- p_base %+%
  filter(drfdata, Function %in% c("Emax", "sigmoid Emax")) +
  geom_text(data = filter(drfdata, Function %in% c("Emax", "sigmoid Emax") & !is.na(pointlabel)),
            size = 10, nudge_y=c(-0.1, -0.15), nudge_x=c(10, 80))

p3 <- p_base %+%
  filter(drfdata) +
  geom_text(data = filter(drfdata, !is.na(pointlabel)),
            size = 10, nudge_y=c(-0.1, -0.15, 0.11), nudge_x=c(10, 80, 2.5))


p1


## -----------------------------------------------------------------------------
#| echo: false
#| warning: false
#| fig.height: 6.5
p2


## -----------------------------------------------------------------------------
#| echo: false
#| warning: false
#| fig.height: 6.5
p3


## ----include=FALSE------------------------------------------------------------
# This is the PATHWAY DRF data by group
pathway <- tibble(
  dose = c(0, 70, 210, 280 * 2),
  group = c(
    "placebo",
    "tezepelumab 70 mg q4w",
    "tezepelumab 210 mg q4w",
    "tezepelumab 280 mg q2w"
  ),
  log_est = log(c(0.67, 0.26, 0.19, 0.22)),
  log_stderr = c(0.10304, 0.17689, 0.22217, 0.19108)
)



## ----echo=FALSE---------------------------------------------------------------
gt(pathway) |>
  fmt_number(columns=c("log_est", "log_stderr"), decimals=3) |>
  opt_stylize(style = 6, color = 'blue')


## ----sigemax_plots------------------------------------------------------------
#| echo: FALSE
#| fig-height: 6
# Plot some example sigEmax curves
sig_curves <- expand_grid(
  Dose = seq(0, 100, length.out = 51),
  Emax = 1,
  E0 = 0,
  ED50 = c(5, 25),
  h = c(1 / 3, 1, 3)
) |>
  mutate(
    `Dose response` = E0 + Emax * Dose^h / (Dose^h + ED50^h),
    text = paste0("ED50=", ED50, "/h=", round(h, 2))
  )

ggplot(
  sig_curves,
  aes(x = Dose, y = `Dose response`, col = text, label = text)
) +
  theme_bw(base_size = 24) +
  geom_line(lwd = 3) +
  scale_y_continuous(labels = scales::percent(seq(0, 1, 0.25))) +
  geom_text_repel(
    data = filter(sig_curves, Dose == 50),
    max.iter = 100,
    size = 6,
    nudge_x = c(35, 30, -16, 30, 10, 5),
    segment.color = NA,
    nudge_y = c(0, -0.02, -0.05, -0.02, -0.03, -0.1)
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")


## ----fit_sigemax--------------------------------------------------------------
#| results: 'hide'
form_sig <- bf(
  log_est | se(log_stderr) ~ E0 + Emax * dose^h / (dose^h + ED50^h),
  nlf(h ~ exp(logh)),
  nlf(ED50 ~ exp(logED50)),
  E0 ~ 1,
  Emax ~ 1,
  logh ~ 1,
  logED50 ~ 1,
  nl = TRUE,
  family = gaussian()
)

prior_sig <- prior(normal(0, 1), nlpar = E0) +
  prior(normal(0, 1), nlpar = logh) +
  prior(normal(0, 1), nlpar = Emax) +
  prior(normal(4, 2), nlpar = logED50)


## ----sigemax_prior------------------------------------------------------------
#| echo: FALSE
#| fig-height: 3
# Plot some example sigEmax curves
expand_grid(
  Dose = seq(0, 100, length.out = 51),
  Emax = 1,
  E0 = 0,
  data.frame(
    ED50 = exp(rnorm(n = 1000, mean = 4, sd = 2)),
    h = exp(rnorm(n = 1000))
  )
) |>
  mutate(
    `Dose response` = E0 + Emax * Dose^h / (Dose^h + ED50^h),
    id = paste0(ED50, "-", h)
  ) |>
  ggplot(aes(x = Dose, y = `Dose response`, group = id)) +
  theme_bw(base_size = 24) +
  geom_line(lwd = 0.5, alpha = 0.05) +
  scale_y_continuous(labels = scales::percent(seq(0, 1, 0.25))) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")


## ----sigemax_prior_wide-------------------------------------------------------
#| echo: FALSE
#| fig-height: 3
# Plot some example sigEmax curves
expand_grid(
  Dose = c(seq(0, 5, 0.1), seq(6, 100, length.out = 51)),
  Emax = 1,
  E0 = 0,
  data.frame(
    ED50 = exp(rnorm(n = 1000, mean = 4, sd = 2)),
    h = exp(rnorm(n = 1000, sd = 3))
  )
) |>
  mutate(
    `Dose response` = E0 + Emax * Dose^h / (Dose^h + ED50^h),
    id = paste0(ED50, "-", h)
  ) |>
  ggplot(aes(x = Dose, y = `Dose response`, group = id)) +
  theme_bw(base_size = 24) +
  geom_line(lwd = 0.5, alpha = 0.05) +
  scale_y_continuous(labels = scales::percent(seq(0, 1, 0.25))) +
  scale_color_brewer(palette = "Dark2") +
  theme(legend.position = "none")


## -----------------------------------------------------------------------------
#| eval: FALSE
# fit_sig <- brm(
#   formula = form_sig,
#   data = pathway,
#   prior = prior_sig,
#   control = list(adapt_delta = 0.99),
#   refresh = 0,
#   seed = 34565
# )


## -----------------------------------------------------------------------------
#| echo: FALSE
#| results: 'hide'
fit_sig <- brm(
  formula = form_sig, 
  data = pathway, 
  prior = prior_sig,
  control = list(adapt_delta = 0.99), 
  seed = 34565,
  refresh = 0
) 


## -----------------------------------------------------------------------------
summary(fit_sig)


## -----------------------------------------------------------------------------
#| echo: false
#| fig.height: 4
compare_prior_posterior <- function(param, fit) {
  prior <- fit$prior |>
    parse_dist(prior) |>
    rename(.variable = nlpar) |>
    filter(.variable == param, source == "user")

  model_param <- paste("b", param, "Intercept", sep = "_")

  post <- tibble(post = as_draws_rvars(fit)[[model_param]])
  ggplot(post, aes(y = 1, xdist = post)) +
    stat_slabinterval() +
    stat_slab(
      data = prior,
      aes(xdist = .dist_obj),
      fill = NA,
      colour = "darkblue"
    ) +
    scale_thickness_shared() +
    ylab(NULL) +
    xlab(param) +
    theme(
      axis.ticks.y = element_blank(),   # removes y-axis ticks
      axis.text.y  = element_blank(),   # removes y-axis text
      axis.title.y = element_blank()    # removes y-axis title
    )
}

cmp_sig <- lapply(prior_sig$nlpar, compare_prior_posterior, fit_sig)

# turns the list of plots into a call of plot1 + plot2 + ... which
# merges the different plots into a single one
Reduce("+", cmp_sig)


## -----------------------------------------------------------------------------
#| eval: false
# tibble(dose = seq(0, 560, length.out = 51), log_stderr = 1) |>
#   tidybayes::add_epred_rvars(object = fit_sig) |>
#   mutate(.delta = .epred - .epred[dose == 0]) |>
#   ggplot(aes(x = dose, ydist = .delta)) +
#   ggdist::stat_lineribbon()


## -----------------------------------------------------------------------------
#| echo: false
#| fig.height: 4
tibble(dose = seq(0, 560, length.out = 51), log_stderr = 1) |>
  add_epred_rvars(object = fit_sig) |>
  mutate(.delta = .epred - .epred[dose == 0]) |>
  ggplot(aes(x = dose, ydist = .delta)) +
  theme_bw(base_size = 24) +
  theme(legend.position = "right") +
  geom_hline(yintercept = 0, color = "darkred", lty = 2) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  geom_point(
    aes(x = dose, y = log_rr),
    inherit.aes = FALSE,
    col = "#E69F00",
    data = tibble(
      dose = c(70, 210, 280 * 2),
      log_rr = log(c(1 - 0.62, 1 - 0.71, 1 - 0.66)),
      log_rr_ucl = log(1 - c(0.42, 0.54, 0.47)),
      log_rr_lcl = log(1 - c(0.75, 0.82, 0.79))
    ),
    size = 3
  ) +
  geom_errorbar(
    aes(x = dose, ymin = log_rr_lcl, ymax = log_rr_ucl),
    inherit.aes = FALSE,
    col = "#E69F00",
    data = tibble(
      dose = c(70, 210, 280 * 2),
      log_rr = log(c(1 - 0.62, 1 - 0.71, 1 - 0.66)),
      log_rr_ucl = log(1 - c(0.42, 0.54, 0.47)),
      log_rr_lcl = log(1 - c(0.75, 0.82, 0.79))
    ),
    size = 1,
    width = 30
  ) +
  ylab("Rate ratio vs. placebo") +
  scale_y_continuous(
    breaks = log(c(1, 0.7, 0.5, 0.35, 0.25, 0.15)),
    labels = c(1, 0.7, 0.5, 0.35, 0.25, 0.15)
  ) +
  scale_x_continuous(breaks = c(0, 70, 210, 280 * 2)) +
  xlab("Total tezepelumab dose [mg/month]")


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\small'
form_mbeta <- bf(
  log_est | se(log_stderr) ~
    E0 +
      Emax *
        (delta1 + delta2)^(delta1 + delta2) /
        (delta1^delta1 * delta2^delta2) *
        (dose / 850)^delta1 *
        (1 - dose / 850)^delta2,
  nlf(delta1 ~ exp(logdelta1)),
  nlf(delta2 ~ exp(logdelta2)),
  E0 ~ 1,
  Emax ~ 1,
  logdelta1 ~ 1,
  logdelta2 ~ 1,
  nl = TRUE,
  family = gaussian()
)

prior_mbeta <- prior(normal(0, 1), nlpar = E0) +
  prior(normal(0, 1), nlpar = Emax) +
  prior(normal(0, 1), nlpar = logdelta1) +
  prior(normal(0, 1), nlpar = logdelta2)



## -----------------------------------------------------------------------------
#| eval: FALSE
# fit_mbeta <- brm(
#   form_mbeta,
#   data = pathway,
#   prior = prior_mbeta,
#   control = list(adapt_delta = 0.99),
#   seed = 24355,
#   refresh = 0
# )


## -----------------------------------------------------------------------------
#| echo: FALSE
#| results: 'hide'
fit_mbeta <- brm(
  form_mbeta, data = pathway, prior = prior_mbeta,
  control = list(adapt_delta = 0.99), seed = 24355, refresh = 0
) 


## -----------------------------------------------------------------------------
#| echo: false
#| fig.height: 8
#| fig.align: center
tibble(dose = seq(0, 560, length.out = 51), log_stderr = 1) |>
  add_epred_rvars(object = fit_mbeta) |>
  mutate(.delta = .epred - .epred[dose == 0]) |>
  ggplot(aes(x = dose, ydist = .delta)) +
  theme_bw(base_size = 24) +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, color = "darkred", lty = 2) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  geom_point(
    aes(x = dose, y = log_rr),
    inherit.aes = FALSE,
    col = "#E69F00",
    data = tibble(
      dose = c(70, 210, 280 * 2),
      log_rr = log(c(1 - 0.62, 1 - 0.71, 1 - 0.66)),
      log_rr_ucl = log(1 - c(0.42, 0.54, 0.47)),
      log_rr_lcl = log(1 - c(0.75, 0.82, 0.79))
    ),
    size = 3
  ) +
  geom_errorbar(
    aes(x = dose, ymin = log_rr_lcl, ymax = log_rr_ucl),
    inherit.aes = FALSE,
    col = "#E69F00",
    data = tibble(
      dose = c(70, 210, 280 * 2),
      log_rr = log(c(1 - 0.62, 1 - 0.71, 1 - 0.66)),
      log_rr_ucl = log(1 - c(0.42, 0.54, 0.47)),
      log_rr_lcl = log(1 - c(0.75, 0.82, 0.79))
    ),
    size = 1,
    width = 30
  ) +
  ylab("Rate ratio compared with placebo") +
  scale_y_continuous(
    breaks = log(c(1, 0.7, 0.5, 0.35, 0.25, 0.15)),
    labels = c(1, 0.7, 0.5, 0.35, 0.25, 0.15)
  ) +
  scale_x_continuous(breaks = c(0, 70, 210, 280 * 2)) +
  xlab("Total tezepelumab dose [mg/month]")


## -----------------------------------------------------------------------------
(loo_mbeta <- loo(fit_mbeta))


## -----------------------------------------------------------------------------
loo_exact_sig <- kfold(fit_sig, folds = "loo", silent=2)
(loo_exact_mbeta <- kfold(fit_mbeta, folds = "loo", silent=2)) 


## -----------------------------------------------------------------------------
loo_compare(loo_exact_sig, loo_exact_mbeta) 


## -----------------------------------------------------------------------------
fit_sig$criteria$loo <- loo_exact_sig
fit_mbeta$criteria$loo <- loo_exact_mbeta
(w_dose <- model_weights(fit_sig, fit_mbeta, weights = "loo")) 


## -----------------------------------------------------------------------------
#| echo: false
dose_df <- data.frame(
  dose = seq(min(pathway$dose), max(pathway$dose), length.out = 101),
  log_stderr = 0
)


## -----------------------------------------------------------------------------
bma_model <- dose_df |>
  add_epred_rvars(fit_sig, value = "SigEmax") |>
  add_epred_rvars(fit_mbeta, value = "ModBeta") |>
  mutate(
    BMA = w_dose["fit_sig"] * SigEmax + w_dose["fit_mbeta"] * ModBeta,
    logRR_SigEmax = SigEmax - SigEmax[dose == 0],
    logRR_ModBeta = ModBeta - ModBeta[dose == 0],
    logRR_BMA = BMA - BMA[dose == 0]
  ) |> select(dose, starts_with("logRR"))
bma_model |>
  gt_preview() |>
  fmt_number(decimals = 2)


## -----------------------------------------------------------------------------
#| echo: false
#| fig.height: 6.5
bma_model |>
  pivot_longer(
    cols = starts_with("logRR"),
    names_to = "Model",
    values_to = ".delta"
  ) |>
  mutate(
    Model = factor(str_remove(Model, "^logRR_"), c("SigEmax", "ModBeta", "BMA"))
  ) |>
  ggplot(aes(x = dose, ydist = .delta)) +
  theme_bw(base_size = 18) +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, color = "darkred", lty = 2) +
  stat_lineribbon() +
  scale_fill_brewer(palette = "Blues") +
  geom_point(
    aes(x = dose, y = log_rr),
    inherit.aes = FALSE,
    col = "#E69F00",
    data = tibble(
      dose = c(70, 210, 280 * 2),
      log_rr = log(c(1 - 0.62, 1 - 0.71, 1 - 0.66)),
      log_rr_ucl = log(1 - c(0.42, 0.54, 0.47)),
      log_rr_lcl = log(1 - c(0.75, 0.82, 0.79))
    ),
    size = 3
  ) +
  geom_errorbar(
    aes(x = dose, ymin = log_rr_lcl, ymax = log_rr_ucl),
    inherit.aes = FALSE,
    col = "#E69F00",
    data = tibble(
      dose = c(70, 210, 280 * 2),
      log_rr = log(c(1 - 0.62, 1 - 0.71, 1 - 0.66)),
      log_rr_ucl = log(1 - c(0.42, 0.54, 0.47)),
      log_rr_lcl = log(1 - c(0.75, 0.82, 0.79))
    ),
    size = 1,
    width = 30
  ) +
  facet_wrap(~Model, ncol = 3, nrow = 1, labeller = label_both) +
  ylab("Rate ratio compared with placebo") +
  scale_y_continuous(
    breaks = log(c(1, 0.7, 0.5, 0.35, 0.25, 0.15)),
    labels = c(1, 0.7, 0.5, 0.35, 0.25, 0.15)
  ) +
  scale_x_continuous(breaks = c(0, 70, 210, 280 * 2)) +
  xlab("Total tezepelumab dose [mg/month]")

