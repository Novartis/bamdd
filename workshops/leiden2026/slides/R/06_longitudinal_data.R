here::i_am("slides/R/06_longitudinal_data.R")

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

## ----echo=FALSE---------------------------------------------------------------
adpasi <- readr::read_csv(here::here("data", "longitudinal.csv"), show_col_types = FALSE) |>
  filter(TRT01P %in% c("PBO", "TRT")) |>
  mutate(AVISIT = factor(AVISIT, paste("Week", c(1, 2 * (1:6)))),
         TRT01P = factor(TRT01P, c("PBO", "TRT")))

pasi_data <- filter(adpasi, PARAMCD == "PASITSCO")


## ----echo = FALSE, fig.width = 10, fig.height = 10/1.62, echo = FALSE---------
adp_with_base <- bind_rows(
  adpasi |>
    filter(PARAMCD == "PASITSCO") |>
    transmute(SUBJID, TRT01P, AVISIT, AVAL, PCHG),
  adpasi |>
    filter(PARAMCD == "PASITSCO") |>
    transmute(SUBJID, TRT01P, AVISIT = "Baseline", AVAL = BASE, PCHG = 0) |>
    distinct()
) |>
  mutate(AVISIT = factor(AVISIT, c("Baseline", levels(adpasi$AVISIT))))

ggplot(adp_with_base,
       aes(x = TRT01P, y = AVAL, fill = TRT01P)) +
  geom_boxplot() +
  facet_grid(. ~ AVISIT) +
  theme_default(base_size=20) +
  labs(y = "PASI score", x = NULL,
       title = "Boxplots of PASI score by treatment group and visit") +
  scale_fill_discrete("Treatment group") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        legend.position = "bottom") 
  


## ----echo=FALSE---------------------------------------------------------------
head(select(pasi_data, CHG, BASE, TRT01P, AVISIT, AVISITN, SUBJID), 10)


## ----results='hide'-----------------------------------------------------------
fit_long_mean <- brm(
  CHG ~ BASE + TRT01P * AVISIT + (1 | SUBJID), 
  data = pasi_data, seed = 2454
)
fit_long_mean <- add_criterion(fit_long_mean, "loo")


## -----------------------------------------------------------------------------
sum(rstan::get_elapsed_time(fit_long_mean$fit))


## ----eval=F-------------------------------------------------------------------
# pp_check(fit_long_mean)


## ----echo=F-------------------------------------------------------------------
pp_check(fit_long_mean) + xlab("Change in PASI score from Baseline") + theme_default(base_size=24)


## -----------------------------------------------------------------------------
loo(fit_long_mean)


## ----eval=FALSE---------------------------------------------------------------
# conditional_effects(fit_long_mean, "AVISIT:TRT01P")

## ----echo=FALSE---------------------------------------------------------------
ce <- conditional_effects(fit_long_mean, "AVISIT:TRT01P")
ce <- plot(ce, plot = FALSE)[[1]]
ce + theme_default(base_size = 24) + theme(legend.position = "bottom")


## ----results='hide'-----------------------------------------------------------
fit_long_lin <- brm(
  CHG ~ BASE + TRT01P * AVISITN + (1 | SUBJID), 
  data = pasi_data, seed = 2454
)
fit_long_lin <- add_criterion(fit_long_lin, "loo")


## -----------------------------------------------------------------------------
sum(rstan::get_elapsed_time(fit_long_lin$fit))


## -----------------------------------------------------------------------------
summary(fit_long_lin)


## ----eval=FALSE---------------------------------------------------------------
# conditional_effects(fit_long_lin, "AVISITN:TRT01P")


## ----echo=FALSE---------------------------------------------------------------
ce <- conditional_effects(fit_long_lin, "AVISITN:TRT01P")
ce <- plot(ce, plot = FALSE)[[1]]
ce + theme_default(base_size = 24) + theme(legend.position = "bottom")


## ----results='hide'-----------------------------------------------------------
fit_long_lin_vs <- brm(
  CHG ~ BASE + TRT01P * AVISITN + (1 + AVISITN | SUBJID),
  data = pasi_data, seed = 2454
)
fit_long_lin_vs <- add_criterion(fit_long_lin_vs, "loo")


## -----------------------------------------------------------------------------
sum(rstan::get_elapsed_time(fit_long_lin_vs$fit))


## -----------------------------------------------------------------------------
summary(fit_long_lin_vs)


## ----eval=FALSE---------------------------------------------------------------
# conditional_effects(fit_long_lin_vs, "AVISITN:TRT01P")


## ----echo=FALSE---------------------------------------------------------------
ce <- conditional_effects(fit_long_lin_vs, "AVISITN:TRT01P")
ce <- plot(ce, plot = FALSE)[[1]]
ce + theme_default(base_size = 24) + theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
loo_compare(fit_long_lin, fit_long_lin_vs)


## ----results='hide'-----------------------------------------------------------
fit_long_mo <- brm(
  CHG ~ BASE + TRT01P * mo(AVISITN) + 
    (1 + mo(AVISITN) | SUBJID),
  data = pasi_data, seed = 24543
)
fit_long_mo <- add_criterion(fit_long_mo, "loo")


## -----------------------------------------------------------------------------
sum(rstan::get_elapsed_time(fit_long_mo$fit))


## ----eval=FALSE---------------------------------------------------------------
# conditional_effects(fit_long_mo, "AVISITN:TRT01P")


## ----echo=FALSE---------------------------------------------------------------
ce <- conditional_effects(fit_long_mo, "AVISITN:TRT01P")
ce <- plot(ce, plot = FALSE)[[1]]
ce + theme_default(base_size = 24) + theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
loo_compare(fit_long_lin_vs, fit_long_mo)


## ----results='hide'-----------------------------------------------------------
fit_long_quad <- brm(
  CHG ~ BASE + TRT01P * (AVISITN + I(AVISITN^2)) + 
    (1 + AVISITN + I(AVISITN^2) | SUBJID),
  data = pasi_data, seed = 245124
)
fit_long_quad <- add_criterion(fit_long_quad, "loo")


## -----------------------------------------------------------------------------
sum(rstan::get_elapsed_time(fit_long_quad$fit))


## ----eval=FALSE---------------------------------------------------------------
# conditional_effects(fit_long_quad, "AVISITN:TRT01P")


## ----echo=FALSE---------------------------------------------------------------
ce <- conditional_effects(fit_long_quad, "AVISITN:TRT01P")
ce <- plot(ce, plot = FALSE)[[1]]
ce + theme_default(base_size = 24) + theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
cor(pasi_data$AVISITN, pasi_data$AVISITN^2)


## -----------------------------------------------------------------------------
pasi_data <- pasi_data |> 
  mutate(AVISITN_c = AVISITN - 6)


## -----------------------------------------------------------------------------
cor(pasi_data$AVISITN_c, pasi_data$AVISITN_c^2)


## ----results='hide'-----------------------------------------------------------
fit_long_quad_c <- brm(
  CHG ~ BASE + TRT01P * (AVISITN_c + I(AVISITN_c^2)) + 
    (1 + AVISITN_c + I(AVISITN_c^2) | SUBJID),
  data = pasi_data, seed = 24547
)
fit_long_quad_c <- add_criterion(fit_long_quad_c, "loo")


## -----------------------------------------------------------------------------
sum(rstan::get_elapsed_time(fit_long_quad_c$fit))


## ----eval=FALSE---------------------------------------------------------------
# conditional_effects(fit_long_quad_c, "AVISITN_c:TRT01P")


## ----echo=FALSE---------------------------------------------------------------
ce <- conditional_effects(fit_long_quad_c, "AVISITN_c:TRT01P")
ce <- plot(ce, plot = FALSE)[[1]]
ce + theme_default(base_size = 24) + theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
loo_compare(fit_long_mo, fit_long_quad_c)


## ----results='hide'-----------------------------------------------------------
fit_long_gp <- brm(
  CHG ~ BASE + gp(AVISITN_c, by = TRT01P, k=9) + 
    (1 + (AVISITN_c + I(AVISITN_c^2)) | SUBJID),
  data = pasi_data, seed = 32454,
  control = list(adapt_delta = 0.95)
)
fit_long_gp <- add_criterion(fit_long_gp, "loo")


## -----------------------------------------------------------------------------
sum(rstan::get_elapsed_time(fit_long_gp$fit))


## -----------------------------------------------------------------------------
loo(fit_long_gp)


## ----eval=FALSE---------------------------------------------------------------
# conditional_effects(fit_long_gp, "AVISITN_c:TRT01P")


## ----echo=FALSE---------------------------------------------------------------
ce <- conditional_effects(fit_long_gp, "AVISITN_c:TRT01P")
ce <- plot(ce, plot = FALSE)[[1]]
ce + theme_default(base_size = 24) + theme(legend.position = "bottom")


## -----------------------------------------------------------------------------
loo_compare(fit_long_quad_c, fit_long_gp)


## ----results='hide'-----------------------------------------------------------
fit_long_gp_ar <- brm(
  CHG ~ BASE + gp(AVISITN_c, by = TRT01P, k=9) + 
    ar(AVISITN_c, gr = SUBJID, p = 1) +
    (1 + (AVISITN_c + I(AVISITN_c^2)) | SUBJID),
  data = pasi_data, seed = 24235,
  iter = 4000, warmup = 1000,
  control = list(adapt_delta = 0.95)
)
fit_long_gp_ar <- add_criterion(fit_long_gp_ar, "loo")


## -----------------------------------------------------------------------------
sum(rstan::get_elapsed_time(fit_long_gp_ar$fit))


## -----------------------------------------------------------------------------
posterior_summary(fit_long_gp_ar, variable = "ar") |>
  round(3)


## ----eval=FALSE---------------------------------------------------------------
# conditional_effects(fit_long_gp_ar, "AVISITN_c:TRT01P")


## ----echo=FALSE---------------------------------------------------------------
ce <- conditional_effects(fit_long_gp_ar, "AVISITN_c:TRT01P")
ce <- plot(ce, plot = FALSE)[[1]]
ce + theme_default(base_size = 24) + theme(legend.position = "bottom")

