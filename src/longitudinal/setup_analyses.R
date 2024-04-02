setup_analyses <- function(adpasi){
  
  formulas <- tribble(
    ~endpoint,   ~formula_name,   ~paramcd,     ~family,                                             ~formula, ~correl, ~longitudinal,
       "PASI",    "cell-means", "PASITSCO",  gaussian(),                         CHG ~ BASE + TRT01P * AVISIT,    TRUE,          TRUE,
       "PASI",        "linear", "PASITSCO",  gaussian(),                        CHG ~ BASE + TRT01P * AVISITN,    TRUE,          TRUE,
       "PASI",     "quadratic", "PASITSCO",  gaussian(), CHG ~ BASE + TRT01P * AVISITN + TRT01P * AVISITN ^ 2,    TRUE,          TRUE,
       "PASI",            "gp", "PASITSCO",  gaussian(),       CHG ~ BASE + TRT01P + gp(AVISITN, by = TRT01P),    TRUE,          TRUE,
       "PASI",        "raneff", "PASITSCO",  gaussian(),          CHG ~ BASE + TRT01P * AVISIT + (1 | SUBJID),   FALSE,          TRUE,
       "PASI", "cross-section", "PASITSCO",  gaussian(),                                  CHG ~ BASE + TRT01P,   FALSE,         FALSE
  )
  
  datasets <- tribble(
      ~paramcd, ~longitudinal, ~missing_approach,                                               ~analysis_data,
    "PASITSCO",          TRUE,            "drop",                filter(adpasi, PARAMCD == "PASITSCO", DROPFL),
    "PASITSCO",         FALSE,            "drop", filter(adpasi, PARAMCD == "PASITSCO", DROPFL, AVISITN == 12)
  )
  
  autocor <- tribble(
    ~correl, ~autocor_name, ~autocor,
       TRUE,         "AR1", ~ ar(time = AVISIT, gr = SUBJID, p = 1),
       TRUE,        "COSY", ~ cosy(time = AVISIT, gr = SUBJID)
  )
  
  analyses <- formulas %>%
    full_join(datasets, c("paramcd", "longitudinal"), multiple = "all") %>%
    full_join(autocor, c("correl"), multiple = "all") %>%
    replace_na(list(autocor_name = "none"))
  
  analyses
  
}


