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
   scale_y_continuous(breaks = c(0,1))


## ----eval = FALSE,  echo = TRUE-----------------------------------------------
## y ~ normal(b0[j] + b1 * x, sigma)
## b0[j] ~ normal(b0, sd0)


## ----warning=FALSE, message=FALSE---------------------------------------------
#| fig-height: 3
#| echo: false
inv_logit <- function(x) 1 / (1 + exp(-x))
x <- runif(10000, -10, 10)
gg1 <- ggplot(data.frame(x), aes(x, after_stat(density))) + 
   geom_histogram() +
   theme(axis.title.y = element_blank())
gg2 <- ggplot(data.frame(x), aes(inv_logit(x), after_stat(density))) + 
   geom_histogram() +
   theme(axis.title.y = element_blank())
gg1 + gg2


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
    gather("Component", "value") %>%
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
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
}


## ----conflict, warning=FALSE--------------------------------------------------
#| echo: false
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
gg_conflict1


## ----conflict2, warning=FALSE-------------------------------------------------
#| echo: false
gg_conflict2


## ----conflict3, warning=FALSE-------------------------------------------------
#| echo: false
gg_conflict3


## ----conflict4, warning=FALSE-------------------------------------------------
#| echo: false
gg_conflict4

