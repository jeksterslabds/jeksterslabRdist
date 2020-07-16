#' ---
#' title: "Test: Normal - PDF"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Normal - PDF}
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
context("Test Normal - PDF.")
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
  "Samnple mean ($\\bar{x}$).",
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
#' ## Probability Density Function
#'
#+ pdf
results_dnorm <- dnorm(
  x = x,
  mean = mu,
  sd = sigma,
  log = FALSE
)
results_normpdf <- normpdf(
  x = x,
  mu = mu,
  sigma = sigma,
  log = FALSE
)
results_dnorm_log <- dnorm(
  x = x,
  mean = mu,
  sd = sigma,
  log = TRUE
)
results_normpdf_log <- normpdf(
  x = x,
  mu = mu,
  sigma = sigma,
  log = TRUE
)
#'
#+ plot_01
plot(
  x = x,
  y = results_normpdf
)
#'
#+ plot_02
plot(
  x = x,
  y = results_dnorm
)
#'
#+ plot_03
plot(
  x = x,
  y = results_normpdf_log
)
#'
#+ plot_04
plot(
  x = x,
  y = results_dnorm_log
)
#'
#' ## Summarize Results
#'
#+ results
knitr::kable(
  x = data.frame(
    x = x,
    dnorm = results_dnorm,
    normpdf = results_normpdf,
    dnorm_log = results_dnorm_log,
    normpdf_log = results_normpdf_log
  ),
  row.names = FALSE
)
#'
#' ## Benchmarking
#'
#+ benchmark
microbenchmark(
  dnorm = dnorm(x = x, mean = mu, sd = sigma, log = FALSE),
  normpdf = normpdf(x = x, mu = mu, sigma = sigma, log = FALSE),
  dnorm_log = dnorm(x = x, mean = mu, sd = sigma, log = TRUE),
  normpdf_log = normpdf(x = x, mu = mu, sigma = sigma, log = TRUE)
)
#'
#' ## testthat
#'
#+ testthat_01, echo=TRUE
test_that("normpdf return the same values as dnorm", {
  expect_equivalent(
    round(
      x = results_dnorm,
      digits = 2
    ),
    round(
      x = results_normpdf,
      digits = 2
    )
  )
})
#'
#+ testthat_02, echo=TRUE
test_that("normpdf log = TRUE return the same values as dnorm log = TRUE", {
  expect_equivalent(
    round(
      x = results_dnorm_log,
      digits = 2
    ),
    round(
      x = results_normpdf_log,
      digits = 2
    )
  )
})
