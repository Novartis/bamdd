here::i_am("/home/beanan1/work/brms_course/slides/R/01_overview.R")
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

## ----echo=FALSE, out.height="2in", out.width="2in"----------------------------
knitr::include_graphics("graphics/brms.png")


## -----------------------------------------------------------------------------
#| eval: false
## formula = y ~ 1 + x + (1 | g)


## -----------------------------------------------------------------------------
#| eval: false
## formula = y ~ 1 + x + (1 + x | g)


## -----------------------------------------------------------------------------
#| eval: false
## formula = y ~ 1 + gp(x) + (1 + x | g)


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
## formula = bf(
##   y ~ 1 + x + (1 | g) + ...,
##   par2 ~ 1 + x + (1 | g) + ...,
##   par3 ~ 1 + x + (1 | g) + ...,
## )


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
## formula = bf(
##   y ~ fun(x, nlpar1, nlpar2),
##   nlpar1 ~ 1 + x + (1 | g) + ...,
##   nlpar2 ~ 1 + (1 | g) + ...,
##   nl = TRUE
## )


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
## family = brmsfamily(
##   family = "<family>", link = "<link>",
##   more_link_arguments
## )


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
## family = brmsfamily(family = "gaussian", link = "identity",
##                     link_sigma = "log")


## -----------------------------------------------------------------------------
#| eval: false
#| echo: true
#| mysize: true
#| size: '\small'
## family = brmsfamily(family = "poisson", link = "log")


## -----------------------------------------------------------------------------
#| mysize: true
#| size: '\tiny'
options(
  # how many processor cores would you like to use?
  mc.cores = 4, 
  # how would you like to access Stan?
  brms.backend = "cmdstanr",
  # cache model binaries
  cmdstanr_write_stan_file_dir=here::here("_brms-cache"),
  # no need to normalize likelihoods
  brms.normalize = FALSE,
  # when you are storing your model to file, 
  # how shall it be updated?
  brms.file_refit = "on_change"  
  # alternatives: "never", "always"
  # use "never" for production
)
# create cache directory if not yet available
dir.create(here::here("_brms-cache"), FALSE)

