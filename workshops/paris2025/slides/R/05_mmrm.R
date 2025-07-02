here::i_am("slides/R/05_mmrm.R")

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
design <- tibble(arm = factor(rep(c("Drug 40 mg (N=50)", "Drug 20 mg (N=50)", "Drug 10 mg (N=50)", "Placebo (N=50)"), each=5),
                   levels=rev(c("Drug 40 mg (N=50)", "Drug 20 mg (N=50)", "Drug 10 mg (N=50)", "Placebo (N=50)"))),
       x=rep(c(1, 3.65, 3.9, 3.65, 1), 4),
       y=rep(c(1, 1, 0.6, 0.2, 0.2), 4)+ rep(3:0, each=5),
       design="Parallel group dose finding")
ggplot(design, aes(x=x, y=y, label=arm, fill=arm)) +
  geom_polygon() +
  theme_void() +
  theme(legend.position="none", panel.border=element_rect(fill = NA),
        strip.background=element_rect(fill="lightgrey"),
        strip.text.x = element_text( margin = margin( b = 5, t = 5) ),
        axis.line.x = element_line(color="black"),
        axis.text.x = element_text(color="black", angle=0, hjust=0.5, size=18),
        axis.ticks.x = element_line(color="black"),
        axis.ticks.length.x = unit(0.2, "cm"),
        axis.title.x = element_blank()) +
  geom_text(data=design |> filter(x==1 & y %in% c(1:10)),
            aes(x=x+0.1, y=y-0.3),
            size=10,
            hjust = 0, vjust=0.75) +
  geom_rect(aes(xmin=0, xmax=0.95, ymin=0.2, ymax=4), fill="darkgrey") +
  geom_text(aes(x=0.45, y=1.9+0.2, label="Run-in"), size=10) +
  scale_fill_manual(values=c("#377eb8", "#fd8d3c", "#f03b20", "#bd0026")) +
  scale_x_continuous(breaks = c(0, 0.975, 1.6, 2.15, 2.9, 3.65),
                     labels = c("     Screening", "Baseline /\nRandomization", "Week 2", "Week 4", "Week 8", "Week 12"))


## -----------------------------------------------------------------------------
#| echo: false
#| fig.height: 5
# Correlation matrix between visits (baseline + 4 post-baseline visits)
corr_matrix <- diag(5)
rho <- c(0.6, 0.48, 0.4, 0.375)
corr_matrix[1,2:5] <- rho[1:4]
corr_matrix[2,3:5] <- rho[1:3]
corr_matrix[3,4:5] <- rho[1:2]
corr_matrix[4,5:5] <- rho[1:1]
corr_matrix[lower.tri(corr_matrix)] <- t(corr_matrix)[lower.tri(corr_matrix)]
 
# Standard deviations by visit (baseline + 4 post-baseline visits)
sds <- sqrt(c(0.75, 0.8, 0.85, 0.95, 1.1))

cov_matrix <- diag(sds) %*% corr_matrix %*% diag(sds)
# print(cov_matrix, digits=2)


## -----------------------------------------------------------------------------
#| echo: false
apply_colnames <- function(x,y){
  colnames(x) <- y
  return(x)
}


## -----------------------------------------------------------------------------
#| echo: false
# Simulate from multivariate normal for control group 
# (before adding treatment effect later)
# We simulate 1000 patients and then apply inclusion criteria and keep the
# first 200 that meet them.
set.seed(4095867)
N <- 1000
Nf <- 200 
simulated_data <- rmvnorm(n=N, mean=rep(0, 5), sigma = cov_matrix) |>
    # turn into tibble
    apply_colnames(c("BASE", paste0("visit", 1:4))) |>
    as_tibble() |>
    # Apply inclusion criteria and keep first 200 patients
    filter(BASE>0) |>
    filter(row_number()<=Nf) |>
    ##filter(row_number()<=200) |>
    # Assign subject ID, treatment group and create missing data
    mutate(USUBJID = row_number(),
           TRT01P = dqsample(x=c(0L, 10L, 20L, 40L), size=Nf, replace=T),
           # Simulate dropouts
           dropp2 = plogis(visit1-2),
           dropp3 = plogis(visit2-2),
           dropp4 = plogis(visit3-2),
           dropv = case_when(runif(n=n())<dropp2 ~ 2L,
                             runif(n=n())<dropp3 ~ 3L,
                             runif(n=n())<dropp4 ~ 4L,
                             TRUE ~ 5L),
           visit2 = ifelse(dropv<=2L, NA_real_, visit2),
           visit3 = ifelse(dropv<=3L, NA_real_, visit3),
           visit4 = ifelse(dropv<=3L, NA_real_, visit4)) |>
    dplyr::select(-dropp2, -dropp3, -dropp4, -dropv) |>
    # Turn data into long-format
    pivot_longer(cols=starts_with("visit"), names_to = "AVISIT", values_to="AVAL") |>
    mutate(
        # Assign visit days
        ADY = case_when(AVISIT=="visit1" ~ 2L*7L,
                        AVISIT=="visit2" ~ 4L*7L,
                        AVISIT=="visit3" ~ 8L*7L,
                        AVISIT=="visit4" ~ 12L*7L),
        # Turn to factor with defined order of visits
        AVISIT = factor(AVISIT, paste0("visit", 1:4)),
        # Assume rising treatment effect over time (half there by week 3) with an 
        # Emax dose response (ED50 = 5 mg)
        AVAL = AVAL + 0.9 * (ADY/7)^3/((ADY/7)^3+3^3) * TRT01P/(TRT01P+5),
        # Change from baseline = value - baseline
        CHG = AVAL - BASE,
        TRT01P=factor(TRT01P)) |>
    relocate(USUBJID, TRT01P, AVISIT, ADY, AVAL, CHG, BASE) |>
    # Discard missing data
    filter(!is.na(AVAL))

simulated_data <- simulated_data |> 
  mutate(USUBJID=factor(USUBJID))


## -----------------------------------------------------------------------------
#| echo: false
#| fig-height: 6
p1 <- simulated_data |>
    bind_rows(simulated_data |>
              group_by(USUBJID) |>
              slice_head(n=1) |>
              mutate(AVISIT="Baseline", ADY=0, AVAL=BASE)
              ) |> arrange(USUBJID, ADY) |>
  ggplot(aes(x=AVISIT, y=AVAL, col=TRT01P)) +
  theme_bw(base_size=24) +
  geom_jitter(height=0, width=0.2) +
  scale_color_manual(values=c("#377eb8", "#fd8d3c", "#f03b20", "#bd0026")) +
  guides(x =  guide_axis(angle = 45)) +
  theme(axis.title.x = element_blank(),
        legend.position="none")

p2 <- simulated_data |>
  ggplot(aes(x=AVISIT, y=CHG, fill=TRT01P)) +
  geom_hline(yintercept=0) +
  geom_jitter(alpha=0.4, size=0.65, col="black", shape=22, 
             position=position_jitterdodge(dodge.width=0.75, 
                                           jitter.width=0.2, jitter.height=0)) +
  geom_boxplot(outlier.alpha = 0, alpha=0.8) +
  theme_bw(base_size=24) +
  scale_fill_manual(values=c("#377eb8", "#fd8d3c", "#f03b20", "#bd0026")) + 
  guides(x =  guide_axis(angle = 45)) +
  theme(legend.position="right",
        axis.title.x = element_blank())

p1 + p2


## -----------------------------------------------------------------------------
#| echo: false
simulated_data |> 
  filter(USUBJID %in% c(3, 9, 13)) |> 
  gt() |> 
  fmt_number(columns=c("AVAL", "CHG", "BASE"), decimals=2) |>
  opt_stylize(style = 6, color = "red", add_row_striping = TRUE) |>
  tab_style(locations = cells_body(rows = c(4, 8, 10)),
            style = cell_borders(sides = "bottom", color = "black", weight = px(1.5),style = "solid"))


## PROC MIXED DATA=simulated_data;
##   CLASS TRT01P AVISIT USUBJID;
##   MODEL CHG ~ TRT01P AVISIT BASE TRT01P*AVISIT AVISIT*BASE
##     / SOLUTION DDFM=KR ALPHA = 0.05;
##   REPEATED AVISIT / TYPE=UN SUBJECT = USUBJID R Rcorr GROUP=TRT01P;
##   LSMEANS TRT01P*AVISIT / DIFFS PDIFF CL OM E;
## RUN;

## -----------------------------------------------------------------------------
#| eval: false
## library(mmrm)
## mmrm_fit <- mmrm(
##   formula = CHG ~ TRT01P + AVISIT + BASE + AVISIT:TRT01P +
##     AVISIT:BASE + us(AVISIT | TRT01P / USUBJID),
##   method = "Kenward-Roger",
##   vcov = "Kenward-Roger-Linear", # to match SAS
##   data = mutate(simulated_data, USUBJID=factor(USUBJID))
## )


## -----------------------------------------------------------------------------
#| eval: true
#| echo: false
library(mmrm)
mmrm_fit <- mmrm(
  formula = CHG ~ TRT01P + AVISIT + BASE + AVISIT:TRT01P + 
    AVISIT:BASE + us(AVISIT | TRT01P / USUBJID),
  method = "Kenward-Roger", 
  data = mutate(simulated_data, USUBJID=factor(USUBJID))
)


## -----------------------------------------------------------------------------
#| echo: FALSE
ex_data <- tibble(treatment = c(rep("PHEN/TPM CR 15/92", 15),
                     rep("PHEN/TPM CR 3.75/23", 15),
                     rep("placebo", 15)),
       week = rep(c(0,4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56), 3),
       chg = c(0, -3.73913043478261, -5.88405797101449, -7.50724637681159,-8.95652173913043, -10.0289855072464, -11.1014492753623, -11.8260869565217, -12.3478260869565, -12.8405797101449, -13.2173913043478, -13.304347826087, -13.3333333333333, -13.4492753623188, -13.1884057971014, 0, -2.57971014492753, -3.50724637681159, -4.17391304347826, -4.98550724637681, -5.56521739130435, -6.14492753623188, -6.55072463768116, -6.72463768115942, -6.89855072463768, -7.18840579710145, -7.15942028985507, -6.98550724637681, -6.84057971014492, -6.7536231884058, 0, -1.07246376811594, -1.47826086956522, -1.88405797101449, -2.08695652173913, -2.31884057971014, -2.43478260869565, -2.60869565217391, -2.69565217391304, -2.78260869565217, -2.92753623188406, -2.66666666666666, -2.46376811594203, -2.34782608695652, 
               -2.20289855072464), 
       stderr = c(NA, 0.594202898550725, 0.347826086956522, 0.434782608695651, 0.521739130434784, 0.594202898550725, 0.666666666666667, 0.710144927536232, 0.826086956521739, 0.869565217391306, 0.898550724637681, 0.942028985507246, 0.956521739130435, 0.971014492753623, 1.02898550724638, NA, 0.449275362318841, 0.449275362318841, 0.579710144927535, 0.652173913043478, 0.797101449275361, 0.956521739130435, 1.02898550724638, 1.13043478260869, 1.15942028985507, 1.23188405797102, 1.31884057971014, 1.3768115942029, 1.44927536231884, 1.42028985507246, NA, 0.333333333333332, 0.289855072463769, 0.391304347826086, 0.492753623188404, 0.623188405797102, 0.724637681159422, 0.782608695652174, 0.840579710144927, 0.869565217391304, 0.956521739130434, 1.01449275362319, 1.05797101449275, 1.10144927536232, 1.10144927536232), 
       lcl = c(NA, -4.90374671545133, -6.56578457433277, -8.35940463095944, -9.97911164410786, -11.1936017879151, -12.4080919317224, -13.2179454382966, -13.9669267698374, -14.5448962184406, -14.97851836292, -15.150690710074, -15.2080814924586, -15.3524287965824, -15.205180331918, NA, -3.46027367421365, -4.3878099060977, -5.31012404900872, -6.26374462759858, -7.12750752390873, -8.01967569535715, -8.56749917249773, -8.94024914194382, -9.17097273569861, -9.60285418385369, -9.74430032743688, -9.68400838451167, -9.68110722397109, -9.53734015195544, NA, -1.72578509629596, -2.04636637233045, -2.65100039974755, -3.05273587644002, -3.54026741065539, -3.85504636560874, -4.14258050964004, -4.34315813193221, -4.48692520394787, -4.80228439100932, -4.65503592634498, -4.53735320103513, -4.5066269974644, -4.36169946123252), 
       ucl = c(NA, -2.57451415411388, -5.20233136769621, -6.65508812266374, -7.93393183415301, -8.86436922657765, -9.79480661900228, -10.4342284747469, -10.7287254040756, -11.1362632018492, -11.4562642457756, -11.4580049420999, -11.4585851742081, -11.5461219280553, -11.1716312622849, NA, -1.69914661564142, -2.62668284752548, -3.0377020379478, -3.70726986515503, -4.00292725869996, -4.27017937710661, -4.53395010286458, -4.50902622037501, -4.62612871357675, -4.7739574103492, -4.57454025227326, -4.28700610824195, -4.00005219631876, -3.96990622485616, NA, -0.419142439935926, -0.910155366799983, -1.11711554228143, -1.12117716703824, -1.09741374876489, -1.01451885178257, -1.07481079470778, -1.04814621589387, -1.07829218735647, -1.05278807275879, -0.678297406988347, -0.390183030848927, -0.189025176448641, -0.0440976402167492))
  ggplot(ex_data, aes(x=week, y=chg, ymin=lcl, ymax=ucl, col=treatment, label=treatment)) +  
  theme_bw(base_size=18) +
  theme(legend.position="none") +
  geom_hline(yintercept=0, alpha=0.5) +
  geom_errorbar(width=0.75) +
  geom_line() +
  geom_point() +
  geom_text(data = ex_data |> filter(week==48), size=6,
            nudge_y=c(-3, -3.25, 3)) +
  xlab("Time since randomization (weeks)") +
  ylab("Percent change in body weight (%)") +
  scale_color_manual(values=c("#e41a1c", "#ff7f00", "#377eb8"))


## -----------------------------------------------------------------------------
contrasts(simulated_data$AVISIT) <- MASS::contr.sdif


## -----------------------------------------------------------------------------
#| echo: true
contrasts(simulated_data$AVISIT) |> MASS::fractions()


## -----------------------------------------------------------------------------
# add the intercept
cmat <- cbind("Intercept" = 1, contrasts(simulated_data$AVISIT))
# compute the inverse matrix
solve(cmat) |> MASS::fractions()


## -----------------------------------------------------------------------------
mmrm_model1 <- bf(
  CHG ~ 1 + AVISIT + BASE + BASE:AVISIT + TRT01P + TRT01P:AVISIT 
    + unstr(time = AVISIT, gr = USUBJID),
  sigma ~ 1 + AVISIT + TRT01P + AVISIT:TRT01P
)


## -----------------------------------------------------------------------------
mmrm_prior1 <- prior(normal(0, 2), class=Intercept) +
    prior(normal(0, 1), class=b) +
    prior(normal(0, log(10.0)/1.64), class=Intercept, dpar=sigma) +
    prior(normal(0, log(2.0)/1.64), class=b, dpar=sigma) +  
    prior(lkj(1), class=cortime)


## -----------------------------------------------------------------------------
#| eval: false
## fit_mmrm1 <- brm(
##   formula = mmrm_model1, data = simulated_data, prior = mmrm_prior1, ...
## )


## -----------------------------------------------------------------------------
#| echo: false
#| include: false
fit_mmrm1 <- brm(
  formula = mmrm_model1,
  prior = mmrm_prior1,
  data = simulated_data,
  seed = 234235
)


## -----------------------------------------------------------------------------
#| eval: false
## emm2 <- emmeans(fit_mmrm1, ~ TRT01P | AVISIT, weights="proportional")


## -----------------------------------------------------------------------------
#| echo: false
emm2 <- emmeans(fit_mmrm1, ~ TRT01P | AVISIT, weights="proportional", 
                at = list(AVISIT = c("visit1", "visit4")))
emm2


## -----------------------------------------------------------------------------
#| eval: false
## emm2 |> as.mcmc() |> summarize_draws()


## -----------------------------------------------------------------------------
contrast(emm2, adjust="none", method="trt.vs.ctrl", ref="TRT01P0")


## -----------------------------------------------------------------------------
#| eval: false
## contrast(emm2, adjust="none", method="trt.vs.ctrl", ref="TRT01P0") |> as.mcmc()


## -----------------------------------------------------------------------------
mmrm_model2 <- bf( 
  CHG ~ 1 + AVISIT + mo(TRT01P) + BASE + mo(TRT01P):AVISIT 
    + BASE:AVISIT + unstr(time = AVISIT, gr = USUBJID),
  sigma ~1 + AVISIT + mo(TRT01P) + mo(TRT01P):AVISIT
)


## -----------------------------------------------------------------------------
#| eval: false
## fit_mmrm2 <- brm(formula = mmrm_model2,
##                  data = simulated_data |> mutate(TRT01P=ordered(TRT01P)), prior = mmrm_prior1, ...)


## -----------------------------------------------------------------------------
#| echo: false
#| include: false
fit_mmrm2 <- brm(
  formula = mmrm_model2,
  data = simulated_data |> mutate(TRT01P=ordered(TRT01P)),
  prior = mmrm_prior1,
  seed = 234235
)


## -----------------------------------------------------------------------------
#| echo: false
map2_dfr(list(mmrm_fit, fit_mmrm1, fit_mmrm2),
         c("Frequentist", "Bayesian", "Bayesian (monotonic)"),
         \(x, y) emmeans(x, ~ TRT01P | AVISIT, weights="proportional") |>
                 contrast(adjust="none", method="trt.vs.ctrl", ref="TRT01P0") |>
           confint() |>
           as_tibble() |>
           mutate(Model=y)) |>
    mutate(Week = c(2L, 4L, 8L, 12L)[as.integer(str_extract(AVISIT, "[0-9]+"))],
         Model = factor(Model, 
                        c("Frequentist", "Bayesian", "Bayesian (monotonic)")),
         Dose = str_remove(str_extract(contrast, "[0-9]+ -"), " -"),
         lower.CL = ifelse(is.na(lower.CL), lower.HPD, lower.CL),
         upper.CL = ifelse(is.na(upper.CL), upper.HPD, upper.CL)) |>
    ggplot(aes(x=Dose, y=estimate, ymin=lower.CL, ymax=upper.CL, col=Model)) +
    geom_hline(yintercept=0, col="darkred", lty=2) +
    theme_bw(base_size=24) +
    theme(legend.position="bottom",
        legend.title = element_text(size=22),
        legend.text = element_text(size=22)) +
    geom_point(position=position_dodge(width=0.4), size=4) +
    geom_errorbar(position=position_dodge(width=0.4), lwd=1) +
    scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9")) +
    facet_wrap(~Week, label=label_both, nrow=1, ncol=4) +
    xlab("Dose [mg]") +
    ylab("Difference to placebo") #+ guides(col=guide_legend(nrow=2, byrow=TRUE))

