# applies patch to brms to avoid unnecessary model binary compilations
patch_brms_recompile <- function() {
    if(packageVersion("brms") <= "2.20.4" & cmdstanr::cmdstan_version() >= "2.29.0") {
        warning("Disabling Stan model canonicalization procedure in brms which forces model recompilations.\nPlease do not use this for production runs!\n")
        # canonicalize Stan model file in accordance with the current Stan version... AND NEVER overwrite anything
        fixed_canonicalize_stan_model <- function(stan_file, overwrite_file = TRUE) {
            cmdstan_mod <- cmdstanr::cmdstan_model(stan_file, compile = FALSE)
            out <- utils::capture.output(
                              cmdstan_mod$format(
                                              canonicalize = list("deprecations", "braces", "parentheses"),
                                              overwrite_file = FALSE, backup = FALSE
                                          )
                          )
            paste0(out, collapse = "\n")
        }
        brms <- asNamespace("brms")
        unlockBinding(".canonicalize_stan_model", brms)
        assignInNamespace(".canonicalize_stan_model", fixed_canonicalize_stan_model, ns="brms", envir=brms)
        assign(".canonicalize_stan_model", fixed_canonicalize_stan_model, envir=brms)
        lockBinding(".canonicalize_stan_model", env=brms)
    }
}
