simulate_fake_data <- function(N = 1E4, G = 1E3, P = 3, seed = 46765875) {
  withr::local_seed(seed)

  # regression coefficients
  beta <- rnorm(P)

  # sampled covariates, group means and fake data
  fake <- matrix(rnorm(N * P), ncol = P)
  dimnames(fake) <- list(NULL, paste0("x", 1:P))

  # fixed effect part and sampled group membership
  fake <- transform(
    as.data.frame(fake),
    theta = fake %*% beta,
    g = sample.int(G, N, replace = TRUE)
  )

  # add random intercept by group
  fake <- merge(fake, data.frame(g = 1:G, eta = rnorm(G)), by = "g")

  # linear predictor
  fake <- transform(fake, mu = theta + eta)

  # sample Poisson data
  fake <- transform(fake, y = rpois(N, exp(mu)))

  # shuffle order of data rows to ensure even distribution of computational effort
  fake <- fake[sample.int(N, N), ]

  # drop not needed row names
  rownames(fake) <- NULL

  return(fake)
}
