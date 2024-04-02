# executes brm with parallelization via clustermq
cmq_brm <- function(..., seed, control=list(adapt_delta=0.9), cores, file, chains=4, .log_worker=FALSE) {
    checkmate::assert_integer(as.integer(seed), lower=1, any.missing=FALSE, len=1)
    brms_global <- options()[grep("^brms", names(options()), value=TRUE)]
    cmdstanr_global <- options()[grep("^cmdstanr", names(options()), value=TRUE)]
    dots <- rlang::enquos(...)
    brms_args <- lapply(dots, rlang::eval_tidy)
    suppressWarnings(model <- do.call(brm, modifyList(brms_args, list(chains=0, cores=1))))
    update_args <- list(object=model, seed=seed, cores=1, chains=1, control=control)
    if(!missing(file)) {
        update_args$file <- sub("\\.rds$", "", file)
    }
    if(!missing(file)) {
        brms_args$file <- sub("\\.rds$", "", file)
    }
    master_lib_paths <- .libPaths()
    update_model <- function(chain_id) {
        .libPaths(master_lib_paths)
        library(brms)
        # in case file is part of extra-arguments, we add here the chain_id
        # to get correct by-chain file caching
        if("file" %in% names(update_args)) {
            update_args <- modifyList(update_args, list(file=paste0(update_args$file, "-", chain_id)))
        }
        if("file" %in% names(brms_args)) {
            brms_args <- modifyList(brms_args, list(file=paste0(brms_args$file, "-", chain_id)))
        }
        update_args$chain_id <- chain_id
        brms_args$chain_id <- chain_id
        # ensure the same brms & cmdstanr global options are set
        options(brms_global)
        options(cmdstanr_global)
        ## for the rstan backend we do an update while for cmdstanr we
        ## have to avoid this for whatever reason
        if(model$backend == "cmdstanr") {
            msg <- capture.output(fit <- do.call(brm, modifyList(brms_args, list(chains=1))))
        } else {
            msg <- capture.output(fit <- do.call(update, update_args))
        }
        list(fit=fit, msg=msg)
    }
    n_jobs <- chains
    backend <- getOption("clustermq.scheduler", "multiprocess")
    if(backend %in% c("multiprocess", "multicore")) {
        n_jobs <- min(chains, getOption("mc.cores", 1))
    }
    cores_per_chain <- 1
    if(!is.null(model$threads$threads)) {
        cores_per_chain <- model$threads$threads
    }
    if(chains == 1 & cores_per_chain == 1) {
        ## looks like a debugging run...avoid clustermq
        return(update_model(1)$fit)
    }
    message("Starting ", chains, " chains with a concurrency of ", n_jobs, " and using ", cores_per_chain, " cores per chain with backend ", backend, "...\n")
    cluster_update <- clustermq::Q(update_model, chain_id=1:chains, n_jobs=n_jobs, export=list(update_args=update_args, brms_args=brms_args, brms_global=brms_global, cmdstanr_global=cmdstanr_global, master_lib_paths=master_lib_paths, model=model),
                                   template=list(cores=cores_per_chain),
                                   log_worker=.log_worker)
    fit <- combine_models(mlist=lapply(cluster_update, "[[", "fit"))
    msg <- lapply(cluster_update, "[[", "msg")
    for(i in seq_len(length(msg))) {
        if(length(msg[[i]]) == 0)
            next
        cat(paste0("Output for chain ", i, ":\n"))
        cat(paste(msg[[i]], collapse="\n"), "\n")
    }
    fit$file <- NULL
    fit
}
