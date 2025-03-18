# Install all required packages

# dependencies <- sort(unique(renv::dependencies()$Package))
# cat(paste0("c(\n", paste(paste0("\"", dependencies, "\""), collapse = ",\n"), "\n)\n"))
# also edit the DESCRIPTION file
# cat(paste(paste0("  ", dependencies), collapse = ",\n"))
c(
"assertthat",
"bayesplot",
"BOIN",
"brms",
"checkmate",
"clustermq",
"cmdstanr",
"coda",
"digest",
"directlabels",
"distributional",
"dplyr",
"dqrng",
"emmeans",
"forcats",
"future",
"GGally",
"ggdist",
"ggfortify",
"ggplot2",
"ggpubr",
"ggrepel",
"ggthemes",
"glue",
"gt",
"here",
"hexbin",
"kableExtra",
"knitr",
"latex2exp",
"lme4",
"lubridate",
"MASS",
"Matrix",
"matrixStats",
"meta",
"mmrm",
"multinma",
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
"rstantools",
"scales",
"simsurv",
"stringr",
"survival",
"survminer",
"tibble",
"tidybayes",
"tidyr",
"tidyverse",
"utils",
"vctrs",
"viridis",
"withr"
)

all_pkgs <- .packages(all.available = TRUE)
missing_dependencies <- dependencies[!dependencies %in% all_pkgs]
if(length(missing_dependencies) > 0) {
  cat("Installing missing dependencies:", missing_dependencies, sep = "\n")
  for(pkg in missing_dependencies) {
    install.packages(pkg)
  }
}
