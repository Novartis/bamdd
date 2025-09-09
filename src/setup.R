## common setup code for each section; enables compilation of each
## section on it's own

knitr::opts_chunk$set(echo = TRUE)

here::i_am("src/setup.R")

# 
    library(future)
    cl <- parallel::makePSOCKcluster(parallelly::availableCores())
    parallel::setDefaultCluster(cl)
    on.exit(parallel::stopCluster(cl))
    # 

## avoid that rstan tries to access ressources on the internet
rstan::rstan_options(javascript = FALSE)

# Set defaults for ggplot2 ----
ggplot2::theme_set( ggplot2::theme_bw(base_size=18) +
                    ggplot2::theme(legend.position = "none"))

scale_colour_discrete <- scale_color_discrete <- function(...) {
  # Alternative: ggsci::scale_color_nejm(...)
  ggplot2::scale_colour_brewer(..., palette="Dark2")
}
scale_fill_discrete <- function(...) {
  # Alternative: ggsci::scale_fill_nejm(...)
  ggplot2::scale_fill_brewer(..., palette="Dark2")
}
scale_colour_continuous <- scale_color_continuous <- function(...) {
  viridis::scale_colour_viridis_c(..., option="turbo")
}
ggplot2::update_geom_defaults("point", list(size=2))
ggplot2::update_geom_defaults("line", list(size=1.5))
# To allow adding label to points e.g. as geom_text_repel(data=. %>% filter(1:n()==n()))
# update_geom_defaults("text_repel", list(label.size = NA, fill = rgb(0,0,0,0),
#                                         segment.color = "transparent", size=6))

# variables we would like to set differently in a CI/CD environment
chains <- as.numeric(Sys.getenv("BRMS_NUM_CHAINS", 4))
iter <- as.numeric(Sys.getenv("BRMS_NUM_ITER", 2000))

brms_cache_dir <- Sys.getenv("BRMS_CACHE_DIR", here::here("_brms-cache"))
dir.create(brms_cache_dir, FALSE)

# use this flag as argument to the cache option of selected knitr
# blocks which take long => gets now controlled from Quarto
# use_knitr_cache <- Sys.getenv("USE_KNITR_CACHE", "yes") == "yes"

##knitr::opts_chunk$set(cache.path = file.path(brms_cache_dir, "knitr-cache-"))


options(brms.iter = iter
       ,brms.chains = chains
       ,brms.backend="cmdstanr"
       ,cmdstanr_write_stan_file_dir=brms_cache_dir
       ,width=120
        )

# Set chapter-specific seed
set.seed(sum(as.integer(digest::digest(knitr::current_input(), raw = T))))

# too old brms versions used with too recent cmdstan versions lead to
# a borken caching of model binaries...this is a brute force fix for
# this

source(here::here("src", "brms_utils.R"))
suppressWarnings(patch_brms_recompile())
