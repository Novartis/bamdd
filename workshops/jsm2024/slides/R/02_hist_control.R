here::i_am("slides/R/02_hist_control.R")
library(here)

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

## ----echo = FALSE, include = FALSE--------------------------------------------
AS_region <- withr::with_seed(63523, bind_cols(AS, region=sample(c("asia", "europe", "north_america"), 8, TRUE)))

## ----echo=FALSE---------------------------------------------------------------
kable(AS_region)


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\small'
form_AS <- bf(r | trials(n) ~ 1 + (1|study), 
              family = binomial("logit"))


## -----------------------------------------------------------------------------
#| eval=FALSE
## get_prior(form_AS, data = AS)


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\small'
bprior_AS <- prior(normal(0, 2), class = "Intercept") +
  prior(normal(0, 1), class = "sd", coef = "Intercept", 
        group = "study")


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\small'
#| results: "hide"
fit_AS <- brm(
  form_AS, data = AS, prior = bprior_AS, 
  seed = 2454
)


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\tiny'
summary(fit_AS)


## -----------------------------------------------------------------------------
AS_new <- data.frame(study = "new_study", n = 1)
pe <- posterior_epred(
  fit_AS, newdata = AS_new, allow_new_levels = TRUE, 
  sample_new_levels = "gaussian"
)
posterior_summary(pe)


## ----fig.height=4, fig.width=7------------------------------------------------
pe_mix <- RBesT::automixfit(pe[, 1], type = "beta")
plot(pe_mix)$mix


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\small'
form_AS_region <- bf(r | trials(n) ~ 1 + (1 | region/study), 
                     family = binomial("logit"))


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\footnotesize'
bprior_AS_region <- prior(normal(0, 2), class="Intercept") +
  prior(normal(0, 0.5), class="sd", coef="Intercept", 
        group="region") +
  prior(normal(0, 0.25), class="sd", coef="Intercept", 
        group="region:study")


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\small'
#| results: 'hide'
fit_AS_region <- brm(
  form_AS_region, data = AS_region, 
  prior = bprior_AS_region, seed = 2341
)


## ----mysize=TRUE, size = "\\tiny"---------------------------------------------
summary(fit_AS_region)


## -----------------------------------------------------------------------------
AS_region_new <- data.frame(study = "new_study_asia", 
                            n = 1, region = "asia")

pe_region <- posterior_epred(
  fit_AS_region, newdata = AS_region_new, 
  allow_new_levels = TRUE, 
  sample_new_levels = "gaussian"
)
posterior_summary(pe_region)


## ----fig.height=4, fig.width=7------------------------------------------------
pe_mix_region <- 
  RBesT::automixfit(pe_region[, 1], type = "beta")
plot(pe_mix_region)$mix

