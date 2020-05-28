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
#'   Values of the random variable \eqn{X}.
#' @param mu Numeric.
#'   Location parameter mean (\eqn{\mu}).
#' @param sigma Numeric.
#'   Positive number.
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

#' Normal - Likelihood
#'
#' Calculates the likelihood of \eqn{X}
#' following a normal distribution.
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
#' @inheritParams normal_pdf
#' @references
#'   [Wikipedia: Normal Distribution](https://en.wikipedia.org/wiki/Normal_distribution)
#'
#'   [Wikipedia: Independent and Identically Distributed](https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables)
#'
#'   [Wikipedia: Likelihood Function](https://en.wikipedia.org/wiki/Likelihood_function)
#' @family normal likelihood functions
#' @export
normal_L <- function(mu,
                     sigma,
                     x) {
  # likelihood
  ((1 / sqrt(2 * pi * sigma^2))^(length(x))) * exp((-1 / (2 * sigma^2)) * sum((x - mu)^2))
}

#' Normal - Log-Likelihood
#'
#' Calculates the log-likelihood of \eqn{X}
#' following a normal distribution.
#'
#' The natural log of the likelihood function for the normal (or Gaussian or Gauss or Laplace–Gauss)
#' distribution is given by
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
#'     \mathcal{l}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right) \\
#'     \mathcal{l}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right)
#'     =
#'     \ln
#'     \left\{
#'       \left(
#'         \frac{1}{\sqrt{2 \pi \sigma^2}}
#'       \right)^n
#'       \exp
#'         \left[
#'           -
#'           \frac{1}{2 \sigma^2}
#'           \sum_{i = 1}^{n}
#'           \left(
#'             x_i - \mu
#'           \right)^2
#'         \right]
#'     \right\} \\
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
#' The negative log-likelihood is given by
#'   \deqn{
#'     -
#'     \mathcal{l}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right)
#'     =
#'     -
#'     \left[
#'       -
#'       \frac{n}{2}
#'       \ln 2 \pi \sigma^2
#'       -
#'       \frac{1}{2 \sigma^2}
#'       \sum_{i = 1}^{n}
#'       \left(
#'         x_i - \mu
#'       \right)^2
#'     \right] \\
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
#' @param neg Logical.
#'   If `TRUE`,
#'   returns,
#'   negative log-likelihood.
#' @inheritParams normal_pdf
#' @inherit normal_L references
#' @family normal likelihood functions
#' @export
normal_ll <- function(mu,
                      sigma,
                      x,
                      neg = TRUE) {
  if (neg) {
    # negative log-likelihood
    return(
      (length(x) / 2) * log(2 * pi * sigma^2) + (1 / (2 * sigma^2)) * sum((x - mu)^2)
    )
  } else {
    # log-likelihood
    return(
      (-length(x) / 2) * log(2 * pi * sigma^2) - (1 / (2 * sigma^2)) * sum((x - mu)^2)
    )
  }
}

#' Normal - Two Log-Likelihood
#'
#' Calculates the two log-likelihood of \eqn{X}
#' following a normal distribution.
#'
#' The two log-likelihood for the normal (or Gaussian or Gauss or Laplace–Gauss)
#' distribution is given by
#'   \deqn{
#'     2
#'     \mathcal{l}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right)
#'     =
#'     2
#'     \left[
#'       -
#'       \frac{n}{2}
#'       \ln 2 \pi \sigma^2
#'       -
#'       \frac{1}{2 \sigma^2}
#'       \sum_{i = 1}^{n}
#'       \left(
#'         x_i - \mu
#'       \right)^2
#'     \right] \\
#'     =
#'     -
#'     n
#'     \ln 2 \pi \sigma^2
#'     -
#'     \frac{1}{\sigma^2}
#'     \sum_{i = 1}^{n}
#'     \left(
#'       x_i - \mu
#'     \right)^2
#'   }
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
#' The negative two log-likelihood is given by
#'   \deqn{
#'     -
#'     2
#'     \mathcal{l}
#'       \left(
#'         \mu,
#'         \sigma^2
#'         \mid
#'         x
#'       \right)
#'     =
#'     -
#'     2
#'     \left[
#'       -
#'       \frac{n}{2}
#'       \ln 2 \pi \sigma^2
#'       -
#'       \frac{1}{2 \sigma^2}
#'       \sum_{i = 1}^{n}
#'       \left(
#'         x_i - \mu
#'       \right)^2
#'     \right] \\
#'     =
#'     n
#'     \ln 2 \pi \sigma^2
#'     +
#'     \frac{1}{\sigma^2}
#'     \sum_{i = 1}^{n}
#'     \left(
#'       x_i - \mu
#'     \right)^2 .
#'   }
#'
#' @param neg Logical.
#'   If `TRUE`,
#'   returns,
#'   negative two log-likelihood.
#' @inheritParams normal_pdf
#' @inherit normal_L references
#' @family normal likelihood functions
#' @export
normal_2ll <- function(mu,
                       sigma,
                       x,
                       neg = TRUE) {
  if (neg) {
    # negative 2 log-likelihood
    return(
      length(x) * log(2 * pi * sigma^2) + ((1 / sigma^2) * sum((x - mu)^2))
    )
  } else {
    # 2 log-likelihood
    return(
      -length(x) * log(2 * pi * sigma^2) - ((1 / sigma^2) * sum((x - mu)^2))
    )
  }
}

#' Normal Distribution - Objective Function
#'
#' Objective function to minimize/maximize
#' to estimate parameters of the normal distribution.
#' See [`normal_2ll()`].
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param theta Vector of parameters \eqn{\theta} of the normal distribution
#'   (`theta[1] = mu` (\eqn{\mu}) and `theta[2] = sigma` (\eqn{\sigma})).
#' @inheritParams normal_ll
#' @inherit normal_L references
#' @family normal likelihood functions
#' @export
normal_obj <- function(theta,
                       x,
                       neg = TRUE) {
  if (theta[2] < 0) {
    return(NA)
  }
  normal_2ll(
    mu = theta[1],
    sigma = theta[2],
    x = x,
    neg = neg
  )
}
