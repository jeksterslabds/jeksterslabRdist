#' Normal Distribution - Probablity Density Function
#'
#' @description Calculates probablities
#' from the probability density function
#' of the normal distribution
#' \eqn{
#'   X
#'   \sim
#'   \mathcal{N}
#'   \left(
#'     \mu,
#'     \sigma^2
#'   \right)
#'   %(\#eq:dist-X-norm)
#' } .
#' This function is identical to [`dnorm()`].
#'
#' @details The normal (or Gaussian or Gauss or Laplace–Gauss)
#' distribution is given by
#' \deqn{
#'   X
#'   \sim
#'   \mathcal{N}
#'   \left(
#'     \mu,
#'     \sigma^2
#'   \right)
#'   %(\#eq:dist-X-norm)
#' }
#' and has
#' the probability density function (PDF)
#' \deqn{
#'   f
#'   \left(
#'     x
#'   \right)
#'   =
#'   \frac{1}{\sigma \sqrt{2 \pi}}
#'   \exp
#'   \left[
#'     -
#'     \frac{1}{2}
#'     \left(
#'       \frac{x - \mu}{\sigma}
#'     \right)^2
#'   \right]
#'   %(\#eq:dist-normpdf-1)
#' }
#' or
#' \deqn{
#'   f
#'   \left(
#'     x
#'   \right)
#'   =
#'   \frac{1}{\sqrt{2 \pi \sigma^2}}
#'   \exp
#'     \left[
#'       -
#'       \frac{
#'       \left(
#'         x - \mu
#'       \right)^2}
#'       {2 \sigma^2}
#'     \right]
#'   %(\#eq:dist-normpdf-2)
#' }
#' with
#' - \eqn{x \in \mathbf{R}},
#' - \eqn{\mu} is the location parameter mean
#'   \eqn{\left( \mu \in \mathbf{R} \right)}, and
#' - \eqn{\sigma^2} is the scale parameter variance
#'   \eqn{\left( \sigma^2 > 0 \right)}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family normal likelihood functions
#' @keywords normal
#' @param x Numeric vector.
#'   Values of the random variable \eqn{X}.
#' @param mu Numeric.
#'   Location parameter mean \eqn{\mu}.
#' @param sigma Numeric.
#'   Positive number.
#'   Scale parameter standard deviation \eqn{\sigma = \sqrt{\sigma^2}}.
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
normpdf <- function(x,
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
#' @description Calculates the likelihood of \eqn{X}
#' following a normal distribution.
#'
#' @details The likelihood function for the normal (or Gaussian or Gauss or Laplace–Gauss)
#' distribution is given by
#' \deqn{
#'   \mathcal{L}
#'     \left(
#'       \mu,
#'       \sigma^2
#'       \mid
#'       x
#'     \right)
#'   =
#'   \prod_{i = 1}^{n}
#'   \left\{
#'     \frac{1}{\sigma \sqrt{2 \pi}}
#'     \exp
#'       \left[
#'         -
#'         \frac{1}{2}
#'         \left(
#'           \frac{x_i - \mu}{\sigma}
#'         \right)^2
#'       \right]
#'   \right\}
#'   %(\#eq:dist-normL-1)
#' }
#'   or
#' \deqn{
#'   \mathcal{L}
#'     \left(
#'       \mu,
#'       \sigma^2
#'       \mid
#'       x
#'     \right)
#'   =
#'   \prod_{i = 1}^{n}
#'   \left\{
#'     \frac{1}{\sqrt{2 \pi \sigma^2}}
#'     \exp
#'       \left[
#'         -
#'         \frac{
#'         \left(
#'           x_i - \mu
#'         \right)^2}
#'         {2 \sigma^2}
#'       \right]
#'   \right\} \\
#'   =
#'   \left(
#'     \frac{1}{\sqrt{2 \pi \sigma^2}}
#'   \right)^n
#'   \exp
#'     \left[
#'       -
#'       \frac{1}{2 \sigma^2}
#'       \sum_{i = 1}^{n}
#'       \left(
#'         x_i - \mu
#'       \right)^2
#'     \right]
#'   %(\#eq:dist-normL-2)
#' }
#' with
#' independent and identically distributed
#' sample data \eqn{x \in \mathbf{R}},
#' \eqn{\mu}
#' is the location parameter mean being estimated
#' \eqn{\left( \mu \in \mathbf{R} \right)},
#' and
#' \eqn{\sigma^2}
#' is the scale parameter variance being estimated
#' \eqn{\left( \sigma^2 > 0 \right)}.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family normal likelihood functions
#' @keywords normal
#' @inheritParams normpdf
#' @references
#'   [Wikipedia: Normal Distribution](https://en.wikipedia.org/wiki/normaldistribution)
#'
#'   [Wikipedia: IID](https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables)
#'
#'   [Wikipedia: Likelihood Function](https://en.wikipedia.org/wiki/Likelihood_function)
#' @export
normL <- function(mu,
                  sigma,
                  x) {
  # likelihood
  ((1 / sqrt(2 * pi * sigma^2))^(length(x))) * exp((-1 / (2 * sigma^2)) * sum((x - mu)^2))
}

#' Normal - Log-Likelihood
#'
#' @description Calculates the log-likelihood of \eqn{X}
#' following a normal distribution.
#'
#' @details The natural log of the likelihood function for the normal (or Gaussian or Gauss or Laplace–Gauss)
#' distribution is given by
#' \deqn{
#'   \ln
#'   \mathcal{L}
#'     \left(
#'       \mu,
#'       \sigma^2
#'       \mid
#'       x
#'     \right)
#'   =
#'   \mathcal{l}
#'     \left(
#'       \mu,
#'       \sigma^2
#'       \mid
#'       x
#'     \right) \\
#'   \mathcal{l}
#'     \left(
#'       \mu,
#'       \sigma^2
#'       \mid
#'       x
#'     \right)
#'   =
#'   \ln
#'   \left\{
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
#'   \right\} \\
#'   =
#'   -
#'   \frac{n}{2}
#'   \ln 2 \pi \sigma^2
#'   -
#'   \frac{1}{2 \sigma^2}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i - \mu
#'   \right)^2 .
#'   %(\#eq:dist-normll)
#' }
#' with
#' independent and identically distributed
#' sample data \eqn{x \in \mathbf{R}},
#' \eqn{\mu}
#' is the location parameter mean being estimated
#' \eqn{\left( \mu \in \mathbf{R} \right)},
#' and
#' \eqn{\sigma^2}
#' is the scale parameter variance being estimated
#' \eqn{\left( \sigma^2 > 0 \right)}.
#'
#' The negative log-likelihood is given by
#' \deqn{
#'   -
#'   \mathcal{l}
#'     \left(
#'       \mu,
#'       \sigma^2
#'       \mid
#'       x
#'     \right)
#'   =
#'   -
#'   \left[
#'     -
#'     \frac{n}{2}
#'     \ln 2 \pi \sigma^2
#'     -
#'     \frac{1}{2 \sigma^2}
#'     \sum_{i = 1}^{n}
#'     \left(
#'       x_i - \mu
#'     \right)^2
#'   \right] \\
#'   =
#'   \frac{n}{2}
#'   \ln 2 \pi \sigma^2
#'   +
#'   \frac{1}{2 \sigma^2}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i - \mu
#'   \right)^2 .
#'   %(\#eq:dist-normnegll)
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family normal likelihood functions
#' @keywords normal
#' @inheritParams normpdf
#' @inherit normL references
#' @param neg Logical.
#'   If `TRUE`,
#'   returns,
#'   negative log-likelihood.
#' @export
normll <- function(mu,
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
#' @description Calculates the two log-likelihood of \eqn{X}
#' following a normal distribution.
#'
#' @details The two log-likelihood for the normal (or Gaussian or Gauss or Laplace–Gauss)
#' distribution is given by
#' \deqn{
#'   2
#'   \mathcal{l}
#'     \left(
#'       \mu,
#'       \sigma^2
#'       \mid
#'       x
#'     \right)
#'   =
#'   2
#'   \left[
#'     -
#'     \frac{n}{2}
#'     \ln 2 \pi \sigma^2
#'     -
#'     \frac{1}{2 \sigma^2}
#'     \sum_{i = 1}^{n}
#'     \left(
#'       x_i - \mu
#'     \right)^2
#'   \right] \\
#'   =
#'   -
#'   n
#'   \ln 2 \pi \sigma^2
#'   -
#'   \frac{1}{\sigma^2}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i - \mu
#'   \right)^2
#'   %(\#eq:dist-norm2ll)
#' }
#' with
#' independent and identically distributed
#' sample data \eqn{x \in \mathbf{R}},
#' \eqn{\mu}
#' is the location parameter mean being estimated
#' (\eqn{\mu \in \mathbf{R}}),
#' and
#' \eqn{\sigma^2}
#' is the scale parameter variance being estimated
#' (\eqn{\sigma^2 > 0}).
#'
#' The negative two log-likelihood is given by
#' \deqn{
#'  -
#'  2
#'  \mathcal{l}
#'    \left(
#'      \mu,
#'      \sigma^2
#'      \mid
#'      x
#'    \right)
#'  =
#'  -
#'  2
#'  \left[
#'    -
#'    \frac{n}{2}
#'    \ln 2 \pi \sigma^2
#'    -
#'    \frac{1}{2 \sigma^2}
#'    \sum_{i = 1}^{n}
#'    \left(
#'      x_i - \mu
#'    \right)^2
#'  \right] \\
#'  =
#'  n
#'  \ln 2 \pi \sigma^2
#'  +
#'  \frac{1}{\sigma^2}
#'  \sum_{i = 1}^{n}
#'  \left(
#'    x_i - \mu
#'  \right)^2 .
#'  %(\#eq:dist-normneg2ll)
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family normal likelihood functions
#' @keywords normal
#' @inheritParams normpdf
#' @inherit normL references
#' @param neg Logical.
#'   If `TRUE`,
#'   returns,
#'   negative two log-likelihood.
#' @export
norm2ll <- function(mu,
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
#' @description Objective function to minimize/maximize
#' to estimate parameters of the normal distribution.
#' See [`norm2ll()`].
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family normal likelihood functions
#' @keywords normal
#' @inheritParams normll
#' @inherit normL references
#' @param theta Vector of parameters \eqn{\theta} of the normal distribution
#'   (`theta[1] = mu` (\eqn{\mu}) and `theta[2] = sigma` (\eqn{\sigma})).
#' @export
normobj <- function(theta,
                    x,
                    neg = TRUE) {
  if (theta[2] < 0) {
    return(NA)
  }
  norm2ll(
    mu = theta[1],
    sigma = theta[2],
    x = x,
    neg = neg
  )
}
