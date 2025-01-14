here::i_am("slides/R/01_overview.R")

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

## data {
##   int<lower=0> N;
##   array[N] int<lower=0,upper=1> y;
## }
## parameters {
##   real<lower=0,upper=1> theta;
## }
## model {
##   theta ~ beta(1,1);
##   y ~ bernoulli(theta);
## }

## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
# formula <- y ~ 1 + x + (1 | g)


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
# formula <- y ~ 1 + x + (1 + x | g)


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
# formula <- y ~ 1 + gp(x) + (1 + x | g)


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
# formula <- y ~ 1 + gp(x, k=9) + (1 + x | g)


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
# formula <- bf(
#   y ~ 1 + x + (1 | g) + ...,
#   par2 ~ 1 + x + (1 | g) + ...,
#   par3 ~ 1 + x + (1 | g) + ...,
# )


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
# formula <- bf(
#   y ~ fun(x, nlpar1, nlpar2),
#   nlpar1 ~ 1 + x + (1 | g) + ...,
#   nlpar2 ~ 1 + (1 | g) + ...,
#   nl = TRUE
# )


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
# family <- brmsfamily(
#     family = "<family>", link = "<link>",
#     more_link_arguments
# )


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
# family <- brmsfamily(family = "gaussian", link = "identity",
#                      link_sigma = "log")


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
# family <- brmsfamily(family = "poisson", link = "log")


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
# library(brms)
# here::i_am("relative_path_in_project_to/rscript.R")
# library(here)
# 
# options(
#   # how many processor cores would you like to use?
#   mc.cores = 4,
#   # how would you like to access Stan?
#   brms.backend = "cmdstanr",
#   # cache model binaries
#   cmdstanr_write_stan_file_dir = here::here("_brms-cache"),
#   # no need to normalize likelihoods
#   brms.normalize = FALSE
# )
# # create cache directory if not yet available
# dir.create(here::here("_brms-cache"), FALSE)

