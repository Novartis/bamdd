library(dplyr)
library(tidyr)
library(ggplot2)
# library(patchwork)
# library(latex2exp)
theme_set(bayesplot::theme_default())

# additional distributions
dbeta2 <- function(x, mu, nu, ...) {
  shape1 <- mu * nu
  shape2 <- (1-mu) * nu
  dbeta(x, shape1, shape2, ...)
}

dbeta2_prime <- function(x, mu, nu) {
  shape1 <- mu * nu
  shape2 <- (1-mu) * nu
  x^(shape1-1) * (1 + x)^(-shape1 - shape2) / beta(shape1, shape2)
}

dbernoulli <- function(x, prob, log = FALSE) {
  dbinom(x, 1, prob, log = log)
}

dcategorical <- function(x, ...) {
  probs <- c(...)
  probs[x]
}

dhyper2 <- function(x, size, k, prob, log = FALSE) {
  dhyper(x, m = size * prob, n = size * (1 - prob), k = k, log = log)
}

plot_dist <- function(dist, bounds, pars, xtype = c("c", "d"), 
                      prefix = c("d", "p", "q"), parnames = NULL, 
                      package = NULL, ...) {
  xtype <- match.arg(xtype)
  prefix <- match.arg(prefix)
  pos <- -1
  if (!is.null(package)) {
    pos <- asNamespace(package)
  }
  dist_fun <- get(paste0(prefix, dist), pos = pos, mode = "function")
  if (xtype == "c") {
    # continuous
    df <- data.frame(x = seq(bounds[1], bounds[2], 0.001))
  } else if (xtype == "d") {
    # discrete
    df <- data.frame(x = bounds[1]:bounds[2])
  }
  if (!is.null(parnames)) {
    parnames <- paste0(parnames, " = ")
  }
  cnames <- rep(NA, length(pars))
  for (i in seq_along(pars)) {
    tmp <- do.call(dist_fun, c(list(df$x), pars[[i]], list(...)))
    cnames[i] <- paste0("$", parnames, pars[[i]], "$", collapse = ", ")
    df[paste0(parnames, pars[[i]], collapse = ", ")] <- tmp
  }
  df <- df %>%
    gather("pars", "dens", -x) %>%
    mutate(pars = factor(pars, unique(pars)))
  
  gg <- ggplot(df, aes(x, dens, color = pars))
  if (xtype == "c") {
    gg <- gg + geom_line(size = 1)
  } else if (xtype == "d") {
    gg <- gg + 
      geom_linerange(aes(ymin=0, ymax=dens), size = 1) +
      geom_line(size = 0.8, linetype = "dotted", alpha = 0.8)
  }
  gg <- gg + 
    scale_color_viridis_d(labels = unname(latex2exp::TeX(cnames))) +
    labs(x = "x", y = "", color = "") + 
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "bottom",
      legend.text = element_text(size = 10)
    )
  if (prefix == "p") {
    gg <- gg +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      theme(axis.ticks.y = element_line(), 
            axis.text.y = element_text(),
            axis.line.y = element_line())
  } else if (prefix == "q") {
    gg <- gg +
      scale_y_continuous() +
      theme(axis.ticks.y = element_line(), 
            axis.text.y = element_text(),
            axis.line.y = element_line())
  }
  gg
}

histogram <- function(x, bins = 30) {
  xname <- deparse(substitute(x))
  xexpr <- eval(parse(text = paste0("expression(", xname, ")")))
  ggplot(data.frame(x), aes(x, y = ..density..)) +
    geom_histogram(fill = "lightgrey", color = "black", bins = bins) +
    xlab(xexpr) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      axis.title.y = element_blank()
    )
}
