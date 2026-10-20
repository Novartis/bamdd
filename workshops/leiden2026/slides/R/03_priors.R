here::i_am("slides/R/03_priors.R")

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

## ----message=FALSE------------------------------------------------------------
#| fig-width: 7
#| fig-height: 3
#| fig-align: center
#| out-width: 75%
#| echo: false
n=1000
#x <- seq(-2.5,2.5, length.out=n)
x=rnorm(n)
p <- brms:::inv_logit(0+1.2*x)
y <- rbinom(n, 1, p)
dat <- data.frame(x, p, y)
ggplot(dat, aes(x, y)) +
   geom_point() +
   geom_smooth(method = "lm") +
   geom_smooth(method = "glm", method.args = list(family = binomial()), 
               color = "green") +
   scale_y_continuous(breaks = c(0,1)) +
   theme_bw(base_size = 24)


## ----bias_variance_target-----------------------------------------------------
#| echo: false
#| fig-width: 6
#| fig-height: 3.6
# illustrative only: each cross is one estimate from a new data set, the bullseye is theta
rings <- tibble(r = seq(2.5, 0.5, by = -0.5),
                fill = c("white", "grey85", "white", "grey85", "firebrick"))
target <- rings |>
  reframe(t = seq(0, 2 * pi, length.out = 200), .by = c(r, fill)) |>
  mutate(x = r * cos(t), y = r * sin(t))

shots <- withr::with_seed(2026, bind_rows(
  tibble(estimator = "atop(hat(theta)[1], 'unbiased, high variance')",
         x = rnorm(25, 0, 1), y = rnorm(25, 0, 1)),
  tibble(estimator = "atop(hat(theta)[2], 'small bias, low variance')",
         x = rnorm(25, 0.6, 0.2), y = rnorm(25, 0.4, 0.2))))

ggplot(shots, aes(x, y)) +
  # group = -r draws the largest ring first so smaller rings stay visible
  geom_polygon(aes(group = -r, fill = fill), data = target, colour = "grey40") +
  geom_point(shape = 4, size = 3, stroke = 1.2, colour = "royalblue4") +
  scale_fill_identity() +
  coord_equal(xlim = c(-2.6, 2.6), ylim = c(-2.6, 2.6)) +
  facet_wrap(~estimator, nrow = 1, labeller = label_parsed) +
  theme_void(base_size = 20)


## ----model_normal_student, include=FALSE--------------------------------------
model_conflict <- stan_model(here::here("stan", "conflict_normal_student-t.stan"))


## ----fit_conflict1, include=FALSE, results='hide'-----------------------------
#| echo: false
sdata_conflict1 <- list(
  df1 = 1000,
  location1 = -5,
  scale1 = 1,
  df2 = 4,
  location2 = 5,
  scale2 = 1
)
fit_conflict1 <- sampling(
  model_conflict, data = sdata_conflict1,
  chains = 1, iter = 20000, warmup = 1000
)


## ----fit_conflict2, include=FALSE, results='hide'-----------------------------
sdata_conflict2 <- list(
  df1 = 4,
  location1 = -5,
  scale1 = 1,
  df2 = 1000,
  location2 = 5,
  scale2 = 1
)
fit_conflict2 <- sampling(
  model_conflict, data = sdata_conflict2,
  chains = 1, iter = 20000, warmup = 1000
)



## ----fit_conflict3, include=FALSE, results='hide'-----------------------------
sdata_conflict3 <- list(
  df1 = 1000,
  location1 = -5,
  scale1 = 1,
  df2 = 1000,
  location2 = 5,
  scale2 = 1
)
fit_conflict3 <- sampling(
  model_conflict, data = sdata_conflict3,
  chains = 1, iter = 20000, warmup = 1000
)


## ----fit_conflict4, include=FALSE, results='hide'-----------------------------
sdata_conflict4 <- list(
  df1 = 4,
  location1 = -5,
  scale1 = 1,
  df2 = 4,
  location2 = 5,
  scale2 = 1
)
fit_conflict4 <- sampling(
  model_conflict, data = sdata_conflict4,
  chains = 1, iter = 20000, warmup = 1000
)


## -----------------------------------------------------------------------------
#| echo: false
plot_conflict <- function(fit, guide = "legend") {
  pars <- 
    as.data.frame(fit) %>%
    select(starts_with("mu")) %>%
    rename(Likelihood = mu1, Prior = mu2, Posterior = mu) %>%
    pivot_longer(everything(), names_to = "Component", values_to = "value") %>%
    mutate(Component = factor(
      Component, levels = c("Likelihood", "Prior", "Posterior")
    ))

  ggplot(pars, aes(value, fill = Component)) +
    geom_density(alpha = 0.7) +
    xlim(c(-10, 10)) +
    ylim(c(0, 0.6)) +
    xlab("Parameter") +
    ylab("") +
    scale_fill_manual(
      values = c("#440154", "#FDE725", "#21908C"),
      guide = guide
    ) +
    #scale_fill_viridis_d(guide = guide) +
    theme_bw(base_size = 24) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.title = element_text(size = 22),
      legend.text = element_text(size = 22)
    )
}


## ----conflict, warning=FALSE--------------------------------------------------
#| echo: false
#| cache: false
gg_conflict1 <- plot_conflict(fit_conflict1) +
  theme(legend.position="bottom")

gg_conflict2 <- plot_conflict(fit_conflict2) +
  theme(legend.position="bottom")

gg_conflict3 <- plot_conflict(fit_conflict3) +
  theme(legend.position="bottom")

gg_conflict4 <- plot_conflict(fit_conflict4) +
  theme(legend.position="bottom")

# gg_conflict1 / gg_conflict2 / gg_conflict3 / gg_conflict4 


## ----conflict1, warning=FALSE-------------------------------------------------
#| echo: false
#| cache: false
gg_conflict1


## ----conflict2, warning=FALSE-------------------------------------------------
#| echo: false
#| cache: false
gg_conflict2


## ----conflict3, warning=FALSE-------------------------------------------------
#| echo: false
#| cache: false
gg_conflict3


## ----conflict4, warning=FALSE-------------------------------------------------
#| echo: false
#| cache: false
gg_conflict4


## ----warning=FALSE, message=FALSE---------------------------------------------
#| fig-height: 4.5
#| echo: false
inv_logit <- function(x) 1 / (1 + exp(-x))
x <- runif(10000, -10, 10)
gg1 <- ggplot(data.frame(x), aes(x, after_stat(density))) + 
   geom_histogram() +
   theme_bw(base_size = 24) +
   theme(axis.title.y = element_blank())
gg2 <- ggplot(data.frame(x), aes(inv_logit(x), after_stat(density))) + 
   geom_histogram() +
   theme_bw(base_size = 24) +
   theme(axis.title.y = element_blank())
gg1 + gg2


## ----echo=FALSE---------------------------------------------------------------
num_draws <- 1E4
priors_cmp <- tibble(prior=c("wide", "brms", "RBesT"),
                     density=c("N(0,10^2)", "S(3, 0, 2.5^2)", "N(0, 2^2)"),
                     sample=rvar(cbind(rnorm(num_draws, 0, 10),
                                      rstudent_t(num_draws, 3, 0, 2.5),
                                      rnorm(num_draws, 0, 2))))
priors_logit <- priors_cmp |>
    ggplot(aes(y=prior, xdist=sample)) +
    stat_slab() +
    xlab("Intercept\nlinear predictor") +
    ylab(NULL) +
    labs(title="Prior logit scale") +
    coord_cartesian(xlim=c(-30, 30)) +
    theme_bw(base_size = 24)

priors_response <- priors_cmp |>
    ggplot(aes(y=prior, xdist=rdo(inv_logit(sample)))) +
    stat_slab() +
    xlab("Intercept\nresponse scale") +
    ylab(NULL) +
    labs(title="Prior response scale") +
    coord_cartesian(xlim=c(0, 1)) +
    theme_bw(base_size = 24)

bayesplot_grid(priors_logit, priors_response, grid_args=list(nrow=1))


## ----echo=FALSE---------------------------------------------------------------

response_category <- function(r) {
  fct_rev(cut(
    pmin(r, 1 - r),
    c(0, 0.01, 0.05, 0.1, 0.5),
    labels = c("extreme\n0–1% & 99–100%", "very_small\n1–5% & 95–99%", "small\n5–10% & 90–95%", "moderate\n10–90%"),
    include.lowest = TRUE
  ))
}

priors_cmp |>
  unnest_rvars() |>
  ggplot(aes(x = prior, fill = response_category(inv_logit(sample)))) +
  geom_bar(position = "fill") +
  xlab("Prior") +
  ylab("Probability") +
  labs(fill = "Category") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  theme_bw(base_size = 24) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22)
  )


