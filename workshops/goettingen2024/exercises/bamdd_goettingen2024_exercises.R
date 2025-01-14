here::i_am("exercises/bamdd_goettingen2024_exercises.R")



## -----------------------------------------------------------------------------
#| eval: false
# cmdstanr::check_cmdstan_toolchain(fix=T)


## -----------------------------------------------------------------------------
#| eval: false
# # install.packages("remotes")
# remotes::install_github("coatless-mac/macrtools")
# macrtools::macos_rtools_install()


## -----------------------------------------------------------------------------
#| eval: false
# # in case packages are missing, please run:
# install.packages(c("here", "ggplot2", "dplyr", "knitr", "brms", "posterior", "tidybayes", "RBesT", "ggrepel", "patchwork", "ggdist", "withr"))


## -----------------------------------------------------------------------------
#| eval: false
# bamdd_zip <- tempfile(fileext=".zip")
# download.file("https://github.com/Novartis/bamdd/archive/refs/heads/main.zip", bamdd_zip)
# ## extracts web-site into the users home
# unzip(bamdd_zip, exdir=normalizePath("~"))
# browseURL(normalizePath(file.path("~", "bamdd-main")))
# # to install all dependencies needed to build the web-site, please run
# source(file.path("~", "bamdd-main", "src", "install_dependencies.R"))


## ----include=FALSE------------------------------------------------------------
here::i_am("exercises/bamdd_goettingen2024_exercises.qmd")


## ----warning = FALSE, message = FALSE-----------------------------------------
#| cache: false
# uncomment below and set to the filename of you R file for the
# exercises
# here::i_am("your-exercise-file.R")
library(here)
library(ggplot2)
library(dplyr)
library(knitr)
library(brms)
library(posterior)
library(tidybayes)
library(bayesplot)
library(RBesT)
library(ggrepel)
library(patchwork)
library(ggdist)

options(
  # how many processor cores would you like to use?
  mc.cores = 4, 
  # how would you like to access Stan?
  brms.backend = "cmdstanr",
  # cache model binaries
  cmdstanr_write_stan_file_dir=here::here("_brms-cache"),
  # no need to normalize likelihoods
  brms.normalize = FALSE
)

dir.create(here::here("_brms-cache"), FALSE)
set.seed(254335)

# in case rstan is used as backend, consider to enable model caching
# by enabling the following line
# rstan_options(auto_write = TRUE)

# Set defaults for ggplot2 ----
theme_set(theme_bw(base_size=18))

scale_colour_discrete <- function(...) {
  # Alternative: ggsci::scale_color_nejm(...)
  # scale_colour_brewer(..., palette="Set1")
  ggthemes::scale_colour_colorblind(...)
}
scale_fill_discrete <- function(...) {
  # Alternative: ggsci::scale_fill_nejm(...)
  #scale_fill_brewer(..., palette="Set1")
  ggthemes::scale_fill_colorblind(...)
}
scale_colour_continuous <- function(...) {
  scale_colour_viridis_c(..., option="turbo")
}
update_geom_defaults("point", list(size=2))
update_geom_defaults("line", list(size=1.5))



## -----------------------------------------------------------------------------
withr::with_seed(654873,
                 AS_region <- bind_cols(RBesT::AS, region=sample(c("asia", "europe", "north_america"), 8, TRUE)))
kable(AS_region)


## -----------------------------------------------------------------------------
model <- bf(r | trials(n) ~ 1 + (1 | study), family=binomial)

model_prior <- prior(normal(0, 2), class=Intercept) +
    prior(normal(0, 1), class=sd, coef=Intercept, group=study)

map_mc_brms  <- brm(model, AS_region, prior = model_prior,
                    seed = 4767, control=list(adapt_delta=0.99),
                    silent = 2, refresh = 0)
map_mc_brms


## -----------------------------------------------------------------------------
#| eval: false
# region_model_fixed <- bf(r | trials(n) ~ 1 + region + (1 | study), family=binomial)


## -----------------------------------------------------------------------------

AS_region_all <- data.frame(region=c("asia", "europe", "north_america", "other")) |>
    mutate(study=paste("new_study", region, sep="_"), r=0, n=6)



## -----------------------------------------------------------------------------
pp_count_AS_region <- posterior_predict(map_mc_brms)

pp_rate_AS_region <- sweep(pp_count_AS_region, 2, AS_region$n, "/")

head(pp_rate_AS_region)

with(AS_region, ppc_intervals(r / n, pp_rate_AS_region))  +
    scale_x_continuous("Study", breaks=1:nrow(AS_region), labels=AS_region$study) +
    ylab("Responder Rate") +
    coord_flip() +
    theme(legend.position="right",
          # suppress vertical grid lines for better readability of intervals
          panel.grid.major.y = element_blank())

pp_AS_region <- AS_region |>
    add_predicted_rvars(map_mc_brms, value="pp_r") |>
    mutate(pp_rate=pp_r/n)

pp_AS_region

with(pp_AS_region, ppc_intervals(r / n, as_draws_matrix(pp_rate)))  +
    scale_x_continuous("Study", breaks=1:nrow(AS_region), labels=AS_region$study) +
    ylab("Responder Rate") +
    coord_flip() +
    theme(legend.position="right",
          # suppress vertical grid lines for better readability of intervals
          panel.grid.major.y = element_blank())



## -----------------------------------------------------------------------------

AS_region_all <- data.frame(region=c("asia", "europe", "north_america", "other")) |>
    mutate(study=paste("new_study", region, sep="_"), r=0, n=6)

## to get brms to include the other factor level in the model, we have
## to add a fake row with region "other" and n=0
AS_region_2 <- mutate(bind_rows(AS_region, mutate(AS_region_all, n=0)[4,]),
                      region=factor(region, levels=c("asia", "europe", "north_america", "other")))

str(AS_region_2)

model_fixed <- bf(r | trials(n) ~ 1 + region + (1|study), family=binomial)

get_prior(model_fixed, AS_region_2)

model_fixed_prior <- prior(normal(0, 2), class=Intercept) +
    prior(normal(0, 1), class=sd, coef=Intercept, group=study) +
    prior(normal(0, log(2)/1.96), class=b)

fixed_mc_brms  <- brm(model_fixed, AS_region_2, prior=model_fixed_prior,
                      seed=4767,
                      silent = 2, refresh = 0, control=list(adapt_delta=0.99))

fixed_mc_brms

post_regions <- posterior_linpred(fixed_mc_brms,
                                  newdata=AS_region_all,
                                  transform=TRUE,
                                  allow_new_levels=TRUE,
                                  sample_new_levels="gaussian")

head(post_regions)

colnames(post_regions) <- AS_region_all$region
map_region <- list()
for(r in AS_region_all$region) {
    map_region[[r]] <- mixfit(post_regions[,r], type="beta", Nc=3, constrain_gt1=TRUE)
}

kable(bind_rows(lapply(map_region, summary), .id="MAP"), digits=3)

sapply(map_region, ess)



## -----------------------------------------------------------------------------
crohn <- RBesT::crohn

crohn_sigma <- 88
crohn$y.se <- crohn_sigma/sqrt(crohn$n)

library(RBesT)
set.seed(1234)
rbest_normal_map_mcmc <- gMAP(cbind(y, y.se) ~ 1 | study, 
                              weights=n, data=crohn,
                              family=gaussian,
                              beta.prior=cbind(0, crohn_sigma),
                              tau.dist="HalfNormal",tau.prior=cbind(0,crohn_sigma/2))

model_normal <- bf(y | se(y.se) ~ 1 + (1 | study), family=gaussian())

prior_normal <- prior(normal(0, 88), class=Intercept) +
    prior(normal(0, 88/2), class=sd, coef=Intercept, group=study)

brms_normal_map_mcmc <- brm(model_normal, crohn, prior=prior_normal,
                            seed=4767,
                            silent = 2, refresh = 0, control=list(adapt_delta=0.99))

## comparing the outputs we see that the random effect posterior
## matches...
rbest_normal_map_mcmc
brms_normal_map_mcmc

brms_normal_map <- posterior_epred(brms_normal_map_mcmc,
                                   newdata=data.frame(study="new", n=1, y=0, y.se=88),
                                   allow_new_levels=TRUE,
                                   sample_new_levels="gaussian")

## ... and the MAP prior is also the same
summarise_draws(brms_normal_map, mean, sd, ~quantile2(., probs = c(0.025, 0.5, 0.975)))



## ----peanutdata---------------------------------------------------------------
#| code-fold: true
#| code-summary: "Show data setup code"
peanut <- tibble(
  TRT01P = structure(c(
    6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L,
    6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 5L, 5L, 5L, 5L, 5L, 5L,
    5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L,
    5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L,
    4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L,
    4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
    3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L,
    3L, 3L, 3L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L
  ),levels = c("PBO", "5 mg", "15 mg", "50 mg", "150 mg", "300 mg"), class = "factor"),
  dose = c(
    300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L,
    300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L,
    300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 300L, 150L, 150L, 150L, 150L, 150L, 150L, 150L,
    150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L,
    150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L, 150L,
    150L, 150L, 150L, 150L, 150L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L,
    50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L,
    50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 50L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L,
    15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L,
    15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 15L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L,
    5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L,
    5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
    0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 0L,
    0L, 0L),
  AVAL = c(
    1000L, 100L, 0L, 600L, 600L, 300L, 300L, 10L, 10L, 100L, 1000L, 300L, 10L, 1000L, 1000L, 1000L, 3L, 3L, 600L, 300L,
    1000L, 1000L, 100L, 1000L, 600L, 600L, 1000L, 1000L, 1000L, 300L, 600L, 1000L, 1000L, 1000L, 1000L, 1000L, 300L,
    1000L, 600L, 1000L, 300L, 1000L, 600L, 600L, 1000L, 1000L, 600L, 1000L, 1000L, 1000L, 1000L, 600L, 1000L, 1000L,
    1000L, 0L, 600L, 300L, 1000L, 1000L, 1000L, 100L, 1000L, 1000L, 600L, 1000L, 1000L, 1000L, 600L, 300L, 600L, 600L,
    100L, 600L, 1000L, 300L, 10L, 3L, 1000L, 300L, 300L, 1000L, 300L, 300L, 10L, 1000L, 1000L, 10L, 10L, 600L, 3L, 10L,
    600L, 600L, 600L, 1000L, 100L, 1000L, 1000L, 10L, 1000L, 3L, 10L, 1000L, 100L, 1000L, 100L, 10L, 300L, 10L, 100L,
    1000L, 0L, 1000L, 100L, 3L, 10L, 100L, 100L, 300L, 1000L, 1000L, 1000L, 10L, 1000L, 3L, 3L, 600L, 600L, 10L, 3L,
    600L, 600L, 300L, 300L, 1000L, 3L, 3L, 1000L, 10L, 1000L, 1000L, 600L, 100L, 300L, 600L, 10L, 100L, 0L, 100L, 3L,
    0L, 10L, 3L, 600L, 300L, 300L, 300L, 600L, 300L, 100L, 3L, 0L, 10L, 600L, 300L, 10L, 300L, 600L, 1000L, 600L, 600L,
    0L, 0L, 600L, 600L, 600L, 0L, 0L, 300L, 100L, 0L, 10L, 300L, 1000L, 300L, 600L, 600L, 300L, 10L, 600L, 100L, 100L,
    300L, 3L, 3L, 300L, 1000L, 10L, 3L, 100L, 3L, 100L, 100L, 300L, 3L, 3L, 600L, 300L, 3L, 3L, 3L, 300L, 3L, 0L, 10L,
    3L, 300L, 10L, 10L, 600L, 0L, 300L, 600L, 0L, 0L, 100L, 100L, 10L, 100L, 10L, 100L, 600L, 0L, 600L, 0L, 10L, 100L,
    0L, 0L, 0L, 0L, 3L, 10L, 3L, 300L, 600L, 0L, 0L, 300L, 10L, 10L, 100L, 300L, 0L, 0L, 0L, 3L, 0L, 0L, 10L, 0L, 300L,
    3L, 0L, 0L, 0L, 0L, 100L, 3L, 0L, 3L, 10L, 10L, 3L, 0L, 3L, 10L, 3L, 0L, 0L, 3L, 100L, 0L, 0L, 300L, 0L, 0L, 0L,
    0L, 0L, 0L, 0L, 3L, 0L, 3L, 0L, 0L, 0L, 0L),
  PARAM = structure(
    c(
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
      1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    levels = "Tolerated amount of peanut protein (mg)", class = "factor"),
  USUBJID = c(
    129L, 103L, 104L, 185L, 34L, 15L, 140L, 151L, 7L, 144L, 200L, 23L, 138L, 115L, 255L, 147L, 224L, 101L, 156L, 284L,
    136L, 248L, 179L, 298L, 168L, 295L, 289L, 241L, 36L, 27L, 123L, 266L, 11L, 291L, 236L, 130L, 173L, 195L, 203L, 160L,
    274L, 167L, 290L, 94L, 9L, 269L, 122L, 135L, 64L, 26L, 95L, 10L, 5L, 234L, 161L, 299L, 88L, 69L, 35L, 233L, 286L,
    85L, 91L, 189L, 80L, 152L, 223L, 287L, 244L, 57L, 108L, 18L, 62L, 157L, 300L, 283L, 164L, 243L, 89L, 220L, 24L, 271L,
    166L, 118L, 201L, 127L, 121L, 41L, 267L, 213L, 49L, 73L, 202L, 134L, 112L, 25L, 227L, 29L, 251L, 273L, 119L, 132L,
    74L, 270L, 83L, 37L, 181L, 258L, 253L, 48L, 120L, 54L, 277L, 176L, 65L, 264L, 107L, 171L, 262L, 162L, 187L, 272L,
    288L, 294L, 245L, 109L, 172L, 204L, 275L, 22L, 66L, 186L, 247L, 17L, 149L, 141L, 177L, 280L, 216L, 40L, 75L, 263L,
    246L, 14L, 81L, 260L, 153L, 45L, 237L, 8L, 117L, 86L, 296L, 146L, 154L, 116L, 38L, 100L, 191L, 175L, 92L, 158L, 192L,
    180L, 256L, 254L, 125L, 222L, 145L, 261L, 155L, 159L, 3L, 206L, 278L, 19L, 31L, 63L, 208L, 55L, 259L, 218L, 111L,
    226L, 33L, 2L, 44L, 297L, 53L, 87L, 16L, 28L, 90L, 207L, 56L, 137L, 128L, 178L, 142L, 143L, 148L, 193L, 229L, 265L,
    97L, 252L, 205L, 150L, 165L, 188L, 52L, 99L, 93L, 1L, 221L, 124L, 210L, 6L, 232L, 21L, 211L, 163L, 96L, 60L, 183L,
    190L, 242L, 42L, 46L, 67L, 126L, 209L, 72L, 194L, 238L, 184L, 39L, 105L, 249L, 61L, 113L, 30L, 77L, 12L, 4L, 51L,
    139L, 20L, 268L, 215L, 292L, 217L, 199L, 32L, 276L, 47L, 225L, 230L, 79L, 71L, 98L, 50L, 13L, 76L, 231L, 250L, 58L,
    68L, 239L, 198L, 293L, 212L, 110L, 59L, 182L, 133L, 170L, 43L, 282L, 281L, 131L, 114L, 196L, 214L, 285L, 70L, 102L,
    279L, 106L, 197L, 257L, 228L, 84L, 235L, 78L, 169L, 240L, 219L, 174L, 82L),
  CRIT1 = structure(c(
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    levels = "Tolerating >=600 mg of peanut protein without dose-limiting symptoms", class = "factor"),
  CRIT1FL = structure(c(
    2L, 1L, 1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L,
    2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L,
    2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L,
    1L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 2L,
    1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 2L, 1L,
    2L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L,
    2L, 2L, 2L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L,
    1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    levels = c("N", "Y"), class = "factor"),
  CRIT2 = structure(c(
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    levels = "Tolerating >=1000 mg of peanut protein without dose-limiting symptoms", class = "factor"),
  CRIT2FL = structure(c(
    2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 2L,
    2L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L, 2L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L,
    1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L,
    1L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 2L, 2L, 1L, 2L, 1L, 1L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 2L,
    1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L,
    2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L),
    levels = c("N", "Y"), class = "factor" ) )


## -----------------------------------------------------------------------------
#| code-fold: true
peanut |>
  ggplot(aes(x=factor(AVAL), fill=TRT01P)) +
  geom_bar( position=position_dodge(preserve = "single")) +
  scale_fill_manual(values=c("royalblue", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d")) +
  xlab("Peanut protein tolerated without dose-limiting symptoms (mg)") +
  ylab("Patients") +
  theme(legend.position="right")

peanut |>
  group_by(TRT01P) |>
  summarize(proportion = sum(CRIT1FL=="Y") / n()) |>
  ggplot(aes(x=TRT01P, y=proportion, fill=TRT01P)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_text(aes(label=scales::percent(proportion)), size=7, vjust=c(0, rep(1.5, 5))) +
  scale_fill_manual(values=c("royalblue", "#fb6a4a", "#ef3b2c", "#cb181d", "#a50f15", "#67000d")) +
  scale_y_continuous(breaks=seq(0, 0.7, 0.1), labels=scales::percent) +
  xlab("Tolerating >=600 mg of peanut protein without dose-limiting symptoms") +
  ylab("Percentage of patients")


## ----peanutmodel1-------------------------------------------------------------
#| eval: false
# model1 <- bf( #---- YOUR CODE HERE --- # ~ E0 + Emax * dose^h/(dose^h + ED50^h),
#               family= #---- YOUR CODE HERE --- # ,
#               nlf(h ~ exp(logh)), nlf(ED50 ~ exp(logED50)),
#               E0 ~ 1, logED50 ~ 1, logh ~ 1, Emax ~ 1,
#               nl=TRUE)
# 
# priors1 <- prior(normal(-2.944, 0.38), nlpar=E0) +
#            prior(normal(0, 1), nlpar=logh) +
#            prior(normal(0, 6), nlpar=Emax) +
#            prior(normal(2.7, 2.5), nlpar=logED50)


## ----fitpeanutmodel1----------------------------------------------------------
#| eval: false
# brmfit1 <- brm(
#   formula = model1,
#   prior = priors1,
#   data = summarize(peanut, y=sum(CRIT1FL=="Y"), n=n(), .by=c(TRT01P, dose))
# )


## -----------------------------------------------------------------------------
#| eval: false
# tibble(dose = seq(0, 300, 1), n=1) |> tidybayes::add_epred_rvars(brmfit1)


## ----sol_peanutmodel1---------------------------------------------------------
#| cache.lazy: false
model1 <- bf( y | trials(n) ~ E0 + Emax * dose^h/(dose^h + ED50^h),
              family=binomial(link="logit"),
              nlf(h ~ exp(logh)), nlf(ED50 ~ exp(logED50)),
              E0 ~ 1, logED50 ~ 1, logh ~ 1, Emax ~ 1, 
              nl=TRUE)

priors1 <- prior(normal(-2.944, 0.38), nlpar=E0) +
           prior(normal(0, 1), nlpar=logh) +
           prior(normal(0, 6), nlpar=Emax) +
           prior(normal(2.7, 2.5), nlpar=logED50)


brmfit1 <- brm(
  formula = model1,
  prior = priors1, 
  data = summarize(peanut, y=sum(CRIT1FL=="Y"), n=n(), .by=c(TRT01P, dose)),
  refresh=0
)

summary(brmfit1)

tibble(dose = seq(0, 300, 1), n=1) |> 
  tidybayes::add_epred_rvars(brmfit1) |>
  ggplot(aes(x=dose, ydist=.epred)) +
  stat_lineribbon() +
  scale_fill_brewer(palette="Blues") +
  ylab("Model predicted proportion") +
  xlab("Dose [mg]")

pr <- tibble(dose = seq(0, 300, 1), n=1) |> 
    tidybayes::add_epred_rvars(brmfit1)

left_join(pr, filter(pr, dose == 0)|> 
              dplyr::select(-dose) |> 
              rename(.pbo=.epred), 
          by="n") |>
    mutate(diff = .epred - .pbo) |>
    ggplot(aes(x=dose, ydist=diff)) +
    stat_lineribbon() +
    scale_fill_brewer(palette="Blues") +
    ylab("Model predicted difference\nin proportion vs. placebo") +
    xlab("Dose [mg]")


## ----sol_lognormal------------------------------------------------------------
#| cache: false
model2 <- bf( logAVAL | cens("interval", logAVALnext) ~ E0 + Emax * dose^h/(dose^h + ED50^h),
              family = gaussian(),
              nlf(h ~ exp(logh)), nlf(ED50 ~ exp(logED50)),
              E0 ~ 1, logED50 ~ 1, logh ~ 1, Emax ~ 1, 
              nl=TRUE )

priors2 <- prior(normal(0, 5), nlpar=E0) + # centered on 1 mg, but
    # very uncertain
    prior(normal(0, 1), nlpar=logh) + # No particular reason
    # to change prior (dose
    # response could be
    # just as steep here as
    # for binomial outcome)
    prior(normal(0, 6.5), nlpar=Emax) + # decent prior probability
    # for >= 1000 mg, still conceivable could reach 300000 No
    # particular reason to change prior here (ED50 could be same as
    # for binomial outcome)
    prior(normal(2.7, 2.5), nlpar=logED50) + # Additional parameter: residual SD, clearly a decent amount of
    # variability, want to
    # allow large SD on
    # the log-scale
    prior(lognormal(5, 10), class=sigma) 

brmfit2 <- brm(
  formula = model2,
  prior = priors2,
  data = peanut |> 
    mutate(logAVAL = ifelse(AVAL==0, -28, log(AVAL)),
           logAVALnext = case_when(AVAL==0 ~ log(3),
                                   AVAL==3 ~ log(10),
                                   AVAL==10 ~ log(100),
                                   AVAL==100 ~ log(300),
                                   AVAL==300 ~ log(600),
                                   AVAL==600 ~ log(1000),
                                   TRUE ~ log(30000))),
  refresh=0
  )

summary(brmfit2)


## ----sol_ordinal--------------------------------------------------------------
#| cache: false
# Note that we omit E0 from the model, because the expected placebo outcome
# is already described by the intercept parameters of the distribution
model3 <- bf(AVAL ~ Emax * dose^h/(dose^h + ED50^h),
             family = cumulative(link = "logit", threshold = "flexible"),
             nlf(h ~ exp(logh)),
             nlf(ED50 ~ exp(logED50)),
             logED50 ~ 1,
             logh ~ 1,
             Emax ~ 1, 
             nl=TRUE)

priors3 <- prior(normal(0, 1), nlpar=logh) + # No particular reason to
                                             # change prior here
           prior(normal(0, 6), nlpar=Emax) + # Prior could be
                                             # different (odds from
                                             # one category to the
                                             # next not necessarily
                                             # the same as for top two
                                             # vs. rest)
           prior(normal(2.7, 2.5), nlpar=logED50) + # No particular
                                                    # reason to change
                                                    # prior here
           prior(student_t(3, 0, 2.5), class=Intercept) +
           prior(student_t(3, -0.2, 1), class=Intercept, coef=1) +
           prior(student_t(3, 2.197, 1), class=Intercept, coef=5)

brmfit3 <- brm(
  formula = model3,
  prior = priors3,
  data = mutate(peanut, AVAL = ordered(AVAL)),
  refresh=0
  )

summary(brmfit3)

