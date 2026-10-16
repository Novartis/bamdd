here::i_am("slides/R/07_time_to_event.R")

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
set.seed(46767)
## n per group
n_grp <- 100
n_hist <- 400
## use month as time-unit
rate_1 <- 1 / 6
rate_cens <- 1 / 10
beta <- c(
  trt = -0.2, ## roughly 20% HR reduction
  soc_alt = 0.05, ## alternative chemotherapy is 5% worse
  hist1 = 0.02
) ## simulated historical data has a 2% worse outcome

## covariates of simulated trial data
covs <- data.frame(
  id = seq(1, 2 * n_grp),
  trt = c(0, 1),
  soc_alt = rbinom(2 * n_grp, 1L, 0.3),
  hist1 = 0L,
  hist2 = 0L
)

## covariates of historical data
hcovs <- data.frame(
  id = seq(2 * n_grp + 1, 2 * n_grp + 1 + n_hist - 1),
  trt = c(0),
  soc_alt = 0,
  hist1 = 1L,
  hist2 = 0L
)

simulate_trial <- function(lambda, gamma, lambda_cens, x, betas) {
  ## simulate censoring times, note that we do not simulate end of
  ## trial with maxt for now...
  cens <- simsurv(lambdas = lambda_cens, gammas = 1, x = x)
  events <- simsurv(lambdas = lambda, gammas = gamma, x = x, betas = betas)
  names(cens) <- paste0(names(cens), "_cens")
  bind_cols(events, select(cens, -id_cens), select(x, -id)) |>
    rename(censtime = eventtime_cens) |>
    mutate(
      event = 1 * (eventtime <= censtime),
      y = if_else(event == 1, eventtime, censtime),
      status = NULL,
      status_cens = NULL
    ) |>
    relocate(id, y, event) |>
    mutate(
      soc = factor(soc_alt, c(0, 1), c("ChA", "ChB")),
      trt_ind = trt,
      trt = factor(trt_ind, c(0, 1), c("ctl", "act")),
      arm = factor(
        paste0(c("ctl", "act")[trt_ind + 1], soc),
        levels = c("actChA", "ctlChA", "actChB", "ctlChB")
      )
    )
}

sim <- simulate_trial(rate_1, 1, rate_cens, covs, beta)
hdata1 <- simulate_trial(rate_1, 1, rate_cens, hcovs, beta) |>
  mutate(id = id + max(sim$id))


## -----------------------------------------------------------------------------
#| echo: false
sim |>
  select(y, event, trt, soc, arm) |>
  gt_preview(8) |>
  gt_format() |>
  fmt_integer(c(event))



## -----------------------------------------------------------------------------
#| eval: false
# cc_inv


## -----------------------------------------------------------------------------
#| echo: false
cc_inv <- matrix(
  c(1/4,  1/4,  1/4,  1/4,
    1/2, -1/2,  1/2, -1/2,
    1/2, -1/2, -1/2,  1/2,
    0,   -1,    0,    1),
  nrow=4, ncol=4, byrow=TRUE,
  dimnames=list(contrast=c("intercept", "effectAvg",
                           "deltaEffect", "deltaControl"),
                arm=c("actChA", "ctlChA",
                      "actChB", "ctlChB")))

cc_inv |> MASS::fractions()


## -----------------------------------------------------------------------------
cc <- solve(cc_inv)


## -----------------------------------------------------------------------------
#| echo: false
cc |> MASS::fractions()


## -----------------------------------------------------------------------------
#| echo: true
contrasts(sim$arm) <- cc[, -1]


## -----------------------------------------------------------------------------
model_weibull1 <- bf(y | cens(1-event) ~ 1 + arm,
                     family=weibull())

prior_weibull1 <- 
  prior(normal(meanInter, log(4)/1.64), class=Intercept) +
  prior(normal(0, sdEffect), coef=armeffectAvg) +
  prior(normal(0, sdDeltaEffect), coef=armdeltaEffect) +
  prior(normal(0, sdDeltaControl), coef=armdeltaControl) +
  prior(gamma(0.1, 0.1), class=shape)


## -----------------------------------------------------------------------------
stanvars_weibull1 <- 
  stanvar(-log(log(2)/8), name = "meanInter") + 
  stanvar(log(2)/1.64, name = "sdEffect") +
  stanvar(log(1.25)/1.64, name = "sdDeltaEffect") +
  stanvar(log(1.25)/1.64, name = "sdDeltaControl")


## -----------------------------------------------------------------------------
#| eval: false
# fit_weibull1 <- brm(
#   formula = model_weibull1,
#   data = sim,
#   prior = weibull_prior1,
#   stanvars = stanvars_weibull1,
#   seed=345665
# )


## -----------------------------------------------------------------------------
#| echo: false
#| include: false
fit_weibull1 <- brm(
  formula = model_weibull1, 
  data = sim, 
  prior = prior_weibull1,
  stanvars = stanvars_weibull1,
  seed=345665
)


## -----------------------------------------------------------------------------
summary(fit_weibull1)


## -----------------------------------------------------------------------------
#| echo: true
#| eval: false
# pp_check(
#   fit_weibull1, type = "km_overlay",
#   ndraws = 100
# )


## -----------------------------------------------------------------------------
#| echo: false
#| fig-height: 4
p_full_fup <- pp_check(
  fit_weibull1,
  type="km_overlay", ndraws=100
) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  xlab("Time [month]") + coord_cartesian(xlim=c(0, 35))

p_full_fup + coord_cartesian(xlim=c(0, 15))


## -----------------------------------------------------------------------------
#| echo: false
hdata1 |>
  select(y, event, arm, hist1) |>
  gt_preview(8) |>
  gt_format() |>
  fmt_integer(c(event, hist1))


## -----------------------------------------------------------------------------
sim_comb <- bind_rows(sim, hdata1)
contrasts(sim_comb$arm) <- cc[, -1]

model_weibull2 <- bf(
  y | cens(1-event) ~ 1 + arm + hist1, 
  family=weibull()
)

prior_weibull2 <- prior_weibull1 +
  prior(student_t(6, 0, sdHist), coef = hist1)

stanvars_weibull2 <- stanvars_weibull1 +
  stanvar(log(1.8)/1.64, name = "sdHist")


## -----------------------------------------------------------------------------
#| eval: false
# fit_weibull2 <- brm(
#   formula = model_weibull2,
#   data = sim_comb,
#   prior = weibull_prior2,
#   stanvars = stanvars_weibull2,
#   seed=4567886
# )


## -----------------------------------------------------------------------------
#| echo: false
#| include: false
fit_weibull2 <- brm(
  formula = model_weibull2, 
  data = sim_comb, 
  prior = prior_weibull2,
  stanvars = stanvars_weibull2,
  seed=4567889
)


## -----------------------------------------------------------------------------
summary(fit_weibull2)

