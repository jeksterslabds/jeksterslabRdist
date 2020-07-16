#' ---
#' title: "Test: Normal - Likelihood"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Normal - Likelihood}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ knitr_options, include=FALSE, cache=FALSE
knitr::opts_chunk$set(
  error = TRUE,
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)
#'
#+ setup
library(testthat)
library(microbenchmark)
library(jeksterslabRdist)
context("Test Normal - Likelihood.")
#'
#' ## Parameters
#'
#+ parameters
n <- 50
mu <- runif(
  n = 1,
  min = 1,
  max = 2
)
sigma <- runif(
  n = 1,
  min = 1,
  max = 2
)
Variable <- c(
  "`n`",
  "`mu`",
  "`sigma`"
)
Description <- c(
  "Sample size ($n$).",
  "Population mean ($\\mu$).",
  "Population standard deviation ($\\sigma$)."
)
Value <- c(
  n,
  mu,
  sigma
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Value
  ),
  row.names = FALSE
)
#'
#' ## Generate Data
#'
#+ generate_data
x <- rnorm(
  n = n,
  mean = mu,
  sd = sigma
)
x_bar <- mean(x)
s <- sd(x)
Variable <- c(
  "`n`",
  "`x_bar`",
  "`s`"
)
Description <- c(
  "Sample size ($n$).",
  "Sample mean ($\\bar{x}$).",
  "Sample standard deviation ($s$)."
)
Value <- c(
  n,
  x_bar,
  s
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Value
  ),
  row.names = FALSE
)
#'
#' ## Likelihood
#'
#+ likelihood
results_dnorm_L <- prod(
  dnorm(
    x = x,
    mean = mu,
    sd = sigma,
    log = FALSE
  )
)
results_dnorm_ll <- -1 * sum(
  dnorm(
    x = x,
    mean = mu,
    sd = sigma,
    log = TRUE
  )
)
results_dnorm_2ll <- -2 * sum(
  dnorm(
    x = x,
    mean = mu,
    sd = sigma,
    log = TRUE
  )
)
results_normL <- normL(
  mu = mu,
  sigma = sigma,
  x = x
)
results_normll <- normll(
  mu = mu,
  sigma = sigma,
  x = x,
  neg = TRUE
)
results_normll_neg_false <- normll(
  mu = mu,
  sigma = sigma,
  x = x,
  neg = FALSE
)
results_normll_L <- -1 * log(results_normL)
results_norm2ll <- norm2ll(
  mu = mu,
  sigma = sigma,
  x = x,
  neg = TRUE
)
results_norm2ll_neg_false <- norm2ll(
  mu = mu,
  sigma = sigma,
  x = x,
  neg = FALSE
)
results_norm2ll_L <- -2 * log(results_normL)
#'
#' ## Maximum Likelihood Estimation
#'
#+ ml
foo <- function(theta, x) {
  -2 * sum(
    dnorm(
      x = x,
      mean = theta[1],
      sd = theta[2],
      log = TRUE
    )
  )
}
results_ml_dnorm <- opt(
  FUN = foo,
  start_values = c(mu, sigma),
  optim = TRUE,
  x = x
)
results_ml_norm2ll <- opt(
  FUN = normobj,
  start_values = c(mu, sigma),
  optim = TRUE,
  x = x,
  neg = TRUE
)
results_nlminb_ml_dnorm <- opt(
  FUN = foo,
  start_values = c(mu, sigma),
  optim = FALSE,
  x = x
)
results_nlminb_ml_norm2ll <- opt(
  FUN = normobj,
  start_values = c(mu, sigma),
  optim = FALSE,
  x = x,
  neg = TRUE
)
#'
#' ## Summarize Results
#'
#+ results
knitr::kable(
  x = data.frame(
    dnorm_ll = results_dnorm_ll,
    normll = results_normll,
    normL = results_normll_L
  ),
  row.names = FALSE,
  caption = "Negative Log-Likelihood"
)
knitr::kable(
  x = data.frame(
    dnorm_2ll = results_dnorm_2ll,
    norm2ll = results_norm2ll,
    normL = results_norm2ll_L
  ),
  row.names = FALSE,
  caption = "Negative Two Log-Likelihood"
)
knitr::kable(
  x = data.frame(
    Item = c("$\\mu$", "$\\sigma$"),
    Parameter = c(mu, sigma),
    Sample = c(x_bar, s),
    dnorm_optim = results_ml_dnorm$par,
    dnorm_nlminb = results_nlminb_ml_dnorm$par,
    normal_negll_optim = results_ml_norm2ll$par,
    normal_negll_nlminb = results_nlminb_ml_norm2ll$par
  ),
  row.names = FALSE,
  col.names = c(
    "Item",
    "Parameter",
    "Closed-form solution",
    "MLE using dnorm (optim)",
    "MLE using dnorm (nlminb)",
    "MLE using normal_negll (optim)",
    "MLE using normal_negll (nlminb)"
  ),
  caption = "Maximum Likelihood Estimates (MLE)"
)
#'
#' ## Benchmarking
#'
#+ benchmark
microbenchmark(
  dnorm_ll = -1 * sum(dnorm(x = x, mean = mu, sd = sigma, log = TRUE)),
  normll = normll(mu = mu, sigma = sigma, x = x, neg = TRUE),
  normll_neg_false = normll(mu = mu, sigma = sigma, x = x, neg = FALSE),
  dnorm_2ll = -2 * sum(dnorm(x = x, mean = mu, sd = sigma, log = TRUE)),
  norm2ll = norm2ll(mu = mu, sigma = sigma, x = x, neg = TRUE),
  norm2ll_neg_false = norm2ll(mu = mu, sigma = sigma, x = x, neg = FALSE)
)
microbenchmark(
  ml_dnorm_optim = opt(FUN = foo, start_values = c(mu, sigma), optim = TRUE, x = x),
  ml_dnorm_nlminb = opt(FUN = foo, start_values = c(mu, sigma), optim = FALSE, x = x),
  ml_normobj_optim = opt(FUN = normobj, start_values = c(mu, sigma), optim = TRUE, x = x, neg = TRUE),
  ml_normobj_nlminb = opt(FUN = normobj, start_values = c(mu, sigma), optim = FALSE, x = x, neg = TRUE)
)
#'
#' ## testthat
#'
#+ testthat_01
test_that("normL returns the same value as dnorm", {
  expect_equivalent(
    round(
      x = results_dnorm_L,
      digits = 2
    ),
    round(
      x = results_normL,
      digits = 2
    )
  )
})
#'
#+ testthat_02
test_that("normll returns the same value as dnorm", {
  expect_equivalent(
    round(
      x = results_dnorm_ll,
      digits = 2
    ),
    round(
      x = results_normll,
      digits = 2
    ),
    round(
      x = results_normll_neg_false * -1,
      digits = 2
    ),
    round(
      x = results_normll_L,
      digits = 2
    )
  )
})
#'
#+ testthat_03
test_that("norm2ll returns the same value as dnorm", {
  expect_equivalent(
    round(
      x = results_dnorm_2ll,
      digits = 2
    ),
    round(
      x = results_norm2ll,
      digits = 2
    ),
    round(
      x = results_norm2ll_neg_false * -1,
      digits = 2
    ),
    round(
      x = results_norm2ll_L,
      digits = 2
    )
  )
})
#'
#+ testthat_04
test_that("normobj maximizes to the same estimates as dnorm", {
  expect_equivalent(
    round(
      x = results_ml_dnorm$par,
      digits = 2
    ),
    round(
      x = results_nlminb_ml_dnorm$par,
      digits = 2
    ),
    round(
      x = results_ml_norm2ll$par,
      digits = 2
    ),
    round(
      x = results_nlminb_ml_norm2ll$par,
      digits = 2
    )
  )
})
#'
#+ testthat_05
test_that("NA output", {
  expect_true(
    is.na(normobj(
      theta = c(-2, -2),
      x = rnorm(10),
      neg = TRUE
    ))
  )
})
