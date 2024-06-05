# Install all required packages

# dependencies <- sort(unique(renv::dependencies()$Package))
dependencies <- c(
  "bayesplot", 
  "BOIN",
  "brms", 
  "checkmate", 
  "clustermq", 
  "cmdstanr", 
  "digest", 
  "directlabels", 
  "dplyr", 
  "dqrng", 
  "emmeans", 
  "future", 
  "GGally", 
  "ggfortify",
  "ggplot2", 
  "ggpubr", 
  "ggrepel", 
  "gt", 
  "here", 
  "knitr", 
  "lme4", 
  "MASS", 
  "Matrix",
  "meta", 
  "mmrm", 
  # "multinma", 
  "mvtnorm", 
  "netmeta",
  "nlme", 
  "OncoBayes2", 
  "parallel", 
  "parallelly", 
  "patchwork", 
  "posterior", 
  "prettyunits",
  "purrr", 
  "RBesT", 
  "readr", 
  "rlang", 
  "rmarkdown", 
  "rstan", 
  "scales", 
  "simsurv", 
  "stringr", 
  "survival", 
  "survminer", 
  "tidybayes", 
  "tidyr", 
  "tidyverse", 
  "utils", 
  "vctrs", 
  "viridis"
)

all_pkgs <- .packages(all.available = TRUE)
missing_dependencies <- dependencies[!dependencies %in% all_pkgs]
if(length(missing_dependencies) > 0) {
  cat("Installing missing dependencies:", missing_dependencies, sep = "\n")
  for(pkg in missing_dependencies) {
    install.packages(pkg)
  }
}
