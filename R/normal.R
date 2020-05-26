#' Normal Distribution - Probablity Density Function
#'
#' Calculates probablities
#' from the probability density function
#' of the normal distribution.
#' This function is identical to [`dnorm()`].
#'
#' The normal (or Gaussian or Gauss or Laplace–Gauss)
#' distribution is given by
#'   \deqn{
#'     \mathcal{N}
#'       \left(
#'         \mu,
#'         \sigma^2
#'       \right)
#'   }
#'   and has
#'   the probability density function (PDF)
#'   \deqn{
#'     f \left( x \right)
#'     =
#'     \frac{1}{\sigma \sqrt{2 \pi}}
#'     \exp
#'       \left[
#'         -
#'         \frac{1}{2}
#'         \left(
#'           \frac{x - \mu}{\sigma}
#'         \right)^2
#'       \right]
#'   }
#'   or
#'   \deqn{
#'     f \left( x \right)
#'     =
#'     \frac{1}{\sqrt{2 \pi \sigma^2}}
#'     \exp
#'       \left[
#'         -
#'         \frac{
#'         \left(
#'           x - \mu
#'         \right)^2}
#'         {2 \sigma^2}
#'       \right]
#'   }
#' with
#' \eqn{x \in \mathbb{R}},
#' \eqn{\mu}
#' is the location parameter mean
#' (\eqn{\mu \in \mathbb{R}}),
#' and
#' \eqn{\sigma^2}
#' is the scale parameter variance
#' (\eqn{\sigma^2 > 0}).
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param x Numeric vector.
#'   Input.
#' @param mu Numeric vector.
#'   Location parameter mean (\eqn{\mu}).
#' @param sigma Numeric vector consisting of positive numbers.
#'   Scale parameter standard deviation (\eqn{\sigma = \sqrt{\sigma^2}}).
#' @param log Logical.
#'   If `TRUE`,
#'   returns \eqn{\log \left( f \left( x \right) \right)}.
#' @return
#'   Returns \eqn{f \left( x \right)}
#'   using the probablity density function
#'   with the supplied parameter/s.
#'   If `log = TRUE`,
#'   returns \eqn{\log \left( f \left( x \right) \right)}.
#' @references
#'   [Wikipedia: Normal Distribution](https://en.wikipedia.org/wiki/Normal_distribution)
#' @export
normal_pdf <- function(x,
                       mu = 0,
                       sigma = 1,
                       log = FALSE) {
  out <- (1 / (sigma * sqrt(2 * pi))) * exp((-1 / 2) * ((x - mu) / sigma)^2)
  if (log) {
    return(log(out))
  }
  out
}

#' Normal - Negative Log-Likelihood
#'
#' Calculates the negative log-likelihood of \eqn{x}
#' following a normal distribution.
#' This function is identical to
#' `-sum(dnorm(x = x, mean = theta[1], sd = theta[2], log = TRUE))`.
#'
#' The likelihood function for the normal (or Gaussian or Gauss or Laplace–Gauss)
#' distribution is given by
#'   \deqn{
#'     \mathcal{L}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right)
#'     =
#'     \prod_{i = 1}^{n}
#'     \left\{
#'       \frac{1}{\sigma \sqrt{2 \pi}}
#'       \exp
#'         \left[
#'           -
#'           \frac{1}{2}
#'           \left(
#'             \frac{x_i - \mu}{\sigma}
#'           \right)^2
#'         \right]
#'     \right\}
#'   }
#'   or
#'   \deqn{
#'     \mathcal{L}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right)
#'     =
#'     \prod_{i = 1}^{n}
#'     \left\{
#'       \frac{1}{\sqrt{2 \pi \sigma^2}}
#'       \exp
#'         \left[
#'           -
#'           \frac{
#'           \left(
#'             x_i - \mu
#'           \right)^2}
#'           {2 \sigma^2}
#'         \right]
#'     \right\} \\
#'     =
#'     \left(
#'       \frac{1}{\sqrt{2 \pi \sigma^2}}
#'     \right)^n
#'     \exp
#'       \left[
#'         -
#'         \frac{1}{2 \sigma^2}
#'         \sum_{i = 1}^{n}
#'         \left(
#'           x_i - \mu
#'         \right)^2
#'       \right]
#'     }
#' with
#' independent and identically distributed
#' sample data \eqn{x \in \mathbb{R}},
#' \eqn{\mu}
#' is the location parameter mean being estimated
#' (\eqn{\mu \in \mathbb{R}}),
#' and
#' \eqn{\sigma^2}
#' is the scale parameter variance being estimated
#' (\eqn{\sigma^2 > 0}).
#'
#' The natural log of the likelihood function is given by
#'   \deqn{
#'     \ln
#'     \mathcal{L}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right)
#'     =
#'     -
#'     \frac{n}{2}
#'     \ln 2 \pi \sigma^2
#'     -
#'     \frac{1}{2 \sigma^2}
#'     \sum_{i = 1}^{n}
#'     \left(
#'       x_i - \mu
#'     \right)^2 .
#'   }
#'
#' The negative log-likelihood is given by
#'   \deqn{
#'     -
#'     \ln
#'     \mathcal{L}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right)
#'     =
#'     \frac{n}{2}
#'     \ln 2 \pi \sigma^2
#'     +
#'     \frac{1}{2 \sigma^2}
#'     \sum_{i = 1}^{n}
#'     \left(
#'       x_i - \mu
#'     \right)^2 .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param theta Vector of parameters \eqn{\theta} of the normal distribution
#'   (`theta[1] = mu` (\eqn{\mu}) and `theta[2] = sigma` (\eqn{\sigma})).
#' @param x Numeric vector.
#'   Independent and identically distributed sample data.
#' @references
#'   [Wikipedia: Normal Distribution](https://en.wikipedia.org/wiki/Normal_distribution)
#'
#'   [Wikipedia: Independent and Identically Distributed Random Variables](https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables)
#'
#'   [Wikipedia: Likelihood Function](https://en.wikipedia.org/wiki/Likelihood_function)
#' @export
normal_negll <- function(theta,
                         x) {
  if (theta[2] < 0) {
    return(NA)
  }
  (length(x) / 2) * log(2 * pi * theta[2]^2) + (1 / (2 * theta[2]^2)) * sum((x - theta[1])^2)
}
