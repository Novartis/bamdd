here::i_am("src/longitudinal/fit_models.R")
library(dplyr)
library(brms)
library(tidyr)
library(ggplot2)
library(readr)
library(here)
library(emmeans)
library(clustermq)
library(purrr)

source(here("src", "longitudinal", "setup_analyses.R"))

# iter <- Sys.getenv("BRMS_ITER")
# chains <- Sys.getenv("BRMS_CHAINS")

adpasi <- readr::read_csv(here("data", "longitudinal.csv")) %>%
  mutate(AVISIT = factor(AVISIT, paste("Week", c(1, 2 * (1:6)))),
         TRT01P = factor(TRT01P, c("PBO", "TRT")))

analyses <- setup_analyses(adpasi)
analyses <- filter(analyses, endpoint == "PASI") %>%
  mutate(id = 1:n())

fit_model <- function(id, analysis_data, formula, formula_name,
                      family, correl, autocor, autocor_name,
                      save_individuals = FALSE){
  
  
  # set up the formula ---------------------------------------------------------
  
  if(correl){
    autocor <- autocor
  } else{
    autocor <- NULL
  }
    
  brms_formula <- bf(
    formula,
    autocor = autocor,
    family = family,
    center = FALSE,
    nl = FALSE
  )
  
  # fit the model --------------------------------------------------------------
  
  prior <- get_prior(
    brms_formula,
    data = analysis_data
  )
  
  cat("\n*** Fitting model", id, ": ", formula_name, ", ", autocor_name, "***\n")
  
  fit <- brm(
    brms_formula,
    prior = prior,
    cores = 4,
    backend = "rstan",
    data = analysis_data
  )
  
  if(save_individuals) saveRDS(fit, here::here("reports", paste0("longitudinal_fit_", id, ".rds")))
  
  # estimate marginal means ----------------------------------------------------
  
  cat("\n*** Estimating marginal means", id, ": ", formula_name, ", ", autocor_name, "***\n")
  
  base <- analysis_data %>%
    select(SUBJID, BASE) %>%
    distinct()
  
  visits <- analysis_data %>%
    select(AVISIT, AVISITN) %>%
    distinct()
  
  treatments <- analysis_data %>%
    select(TRT01P) %>%
    distinct()
  
  emm <- bind_rows(lapply(
    split(visits, 1:nrow(visits)),
    function(visit){
      
      nd <- visit %>%
        crossing(treatments) %>%
        crossing(base) %>%
        split(.$TRT01P)
      
      lp <- map(nd, ~ posterior_linpred(fit, transform = TRUE, newdata = .))
      
      marginal_mean <- lapply(lp, rowMeans) %>%
        as_tibble() %>%
        mutate(diff = .[[2]] - .[[1]]) %>%
        setNames(c(names(nd),
                   paste(rev(names(nd)), collapse = " - ")))
      
      bind_cols(
        visit,
        summarise_draws(as_draws_df(marginal_mean))
      ) %>%
        mutate(
          TRT01P = factor(variable, levels(analysis_data$TRT01P)),
          contrast = case_when(
            is.na(TRT01P) ~ variable,
            TRUE ~ NA_character_
          )
        )
      
    }
  ))
  
  if(save_individuals) saveRDS(emm, here::here("reports", paste0("longitudinal_emm", id, ".rds")))
  
  return(emm)
  
}

emm <- clustermq::Q_rows(
  select(analyses, id, analysis_data, formula, formula_name, family, correl, autocor, autocor_name),
  fit_model,
  const = list(save_individuals = FALSE),
  n_jobs = nrow(analyses),
  pkgs = c("brms", "emmeans", "tidyr", "dplyr", "purrr", "posterior"),
  template = list(
    walltime = 120,
    job_name = "longitudinal_child",
    log_file = here("reports", "longitudinal_child_%I.log"),
    memory = 8000,
    cores = 4
  ),
  job_size = 1
)

saveRDS(emm, file = here("reports", "longitudinal_fits.rds"))
saveRDS(analyses, file = here("reports", "longitudinal_analyses.rds"))

analyses$emm <- emm

