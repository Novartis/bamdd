here::i_am("slides/R/02_hist_control.R")

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

## ----include=FALSE------------------------------------------------------------
#| cache: true
library(RBesT)

gt_format <- function(x, decimals = 2) {
    x |> fmt_number(decimals=decimals)
    ##|>
    ##opt_interactive(
    ##  page_size_default = 6, 
    ##  use_text_wrapping = FALSE, 
    ##  use_compact_mode = TRUE
    ##)
}

map_mc <- gMAP(cbind(r, n-r) ~ 1 | study, family=binomial, data=AS, tau.dist="HalfNormal", tau.prior=1, beta.prior=2)
map_automix <- automixfit(map_mc) 
map_robust  <- robustify(map_automix, mean=0.5, weight=0.2)
pl_AS <- forest_plot(map_mc, size=1.5) + xlab(NULL) + ylab("Responder rate") 
pl_map <- pl_map_AS <- plot(map_mc, size=1.5)$forest_model + xlab(NULL) + ylab("Responder rate")


## flowchart LR
##     M((β,τ)) --> A
##     M --> B
##     M --> C
##     M -->|prediction| T
##     A((θ<sub>1</sub>)) --> YA[Y<sub>1</sub>]
##     B((θ<sub>2</sub>)) --> YB[Y<sub>2</sub>]
##     C((θ<sub>3</sub>)) --> YC[Y<sub>3</sub>]
##     T((θ<sub>*</sub>)) -->|new trial| YT[Y<sub>*</sub>]

## ----fig.width=5--------------------------------------------------------------
#| echo: false
print(pl_AS)


## -----------------------------------------------------------------------------
#| echo: false
summary(map_mc)$theta.pred |> as_tibble() |> gt() |> gt_format()


## ----echo=TRUE,eval=FALSE-----------------------------------------------------
# library(RBesT)
# set.seed(98721487)
# map_mc <- gMAP(cbind(r, n-r) ~ 1 | study,
#                data=RBesT::AS,
#                family=binomial,
#                tau.dist="HalfNormal",
#                tau.prior=1,
#                beta.prior=2)


## -----------------------------------------------------------------------------
#| echo: true
#| eval: true
#| output: false
library(brms)
# use of an outcome modifier to specify number of trials
# intercept study random effect requested with (1 | study)
bmap_model <- bf(r | trials(n) ~ 1 + (1 | study),
                 family=binomial, center=FALSE)
bmap_model_prior <- prior(normal(0, 2), class=b, coef=Intercept) +
    prior(normal(0, 1), class=sd, coef=Intercept, group=study)
bmap_mc <- brm(bmap_model,
               data=RBesT::AS,
               prior=bmap_model_prior,
               seed=98721487,
               control=list(adapt_delta=0.99)) # Stan sampling parameter


## ----echo=FALSE---------------------------------------------------------------
knitr::kable(RBesT::AS)


## -----------------------------------------------------------------------------
#| output: false
form_AS <- bf(r | trials(n) ~ 1 + (1|study), 
              family = binomial("logit"))

get_prior(form_AS, data = AS)

bprior_AS <- prior(normal(0, 2), class = "Intercept") +
    prior(normal(0, 1), class = "sd", coef = "Intercept", group = "study")

fit_AS <- brm(
    form_AS, data = AS, prior = bprior_AS, seed = 2454,
    control=list(adapt_delta=0.99),
    refresh=0
  )


## -----------------------------------------------------------------------------
summary(bmap_mc)


## -----------------------------------------------------------------------------
AS_new <- data.frame(study = "new_study", n = 1)
pe <- posterior_epred(
    bmap_mc, newdata = AS_new,
    allow_new_levels = TRUE, 
    sample_new_levels = "gaussian"
)
posterior_summary(pe)


## -----------------------------------------------------------------------------
pe_mix <- automixfit(pe[, 1], type = "beta")
print(pe_mix, digits=3)


## -----------------------------------------------------------------------------
plot(pe_mix)$mix


## -----------------------------------------------------------------------------
withr::with_seed(654873,
                 AS_region <- bind_cols(AS, region=sample(c("asia", "europe", "north_america"), 8, TRUE)))
kable(AS_region)


## -----------------------------------------------------------------------------
#| output: false
form_AS_region <- bf(r | trials(n) ~ 1 + (1 | region/study), 
                     family = binomial, center=FALSE)

# show which priors brms defines for this model
get_prior(form_AS_region, AS_region)

bprior_AS_region <- prior(normal(0, 2), class=b, coef=Intercept) +
  prior(normal(0, 0.50), class=sd, coef=Intercept, group=region) +
  prior(normal(0, 0.25), class=sd, coef=Intercept, group=region:study)

fit_AS_region <- brm(
  form_AS_region, data = AS_region, prior = bprior_AS_region, seed = 29856341,
  control=list(adapt_delta=0.99))


## -----------------------------------------------------------------------------
summary(fit_AS_region)


## -----------------------------------------------------------------------------
AS_region_new <- data.frame(study = "new_study_asia", 
                            n = 1, region = "asia")

pe_region <- posterior_epred(
  fit_AS_region, newdata = AS_region_new, 
  allow_new_levels = TRUE, 
  sample_new_levels = "gaussian"
)
dim(pe_region)
posterior_summary(pe_region)


## -----------------------------------------------------------------------------
pe_mix_region <- RBesT::automixfit(pe_region[, 1], type = "beta")
plot(pe_mix_region)$mix

