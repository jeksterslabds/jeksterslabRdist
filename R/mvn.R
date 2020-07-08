#' Multivariate Normal - Probablity Density Function
#'
#' Calculates probablities
#' from the probability density function
#' of the multivariate normal distribution
#' \eqn{\mathbf{X} \sim \mathcal{N}_{k} \left( \boldsymbol{\mu}, \boldsymbol{\Sigma}\right)} .
#'
#' The multivariate normal
#' (or multivariate Gaussian,
#' or joint normal)
#' distribution is given by
#' \deqn{
#'   f_{\mathbf{X}}
#'   \left(
#'     X_1,
#'     \dots,
#'     X_k
#'   \right)
#'   =
#'   \frac{
#'     \exp
#'     \left[
#'       -
#'       \frac{1}{2}
#'       \left(
#'         \mathbf{X} - \boldsymbol{\mu}
#'       \right)^{\prime}
#'       \boldsymbol{\Sigma}^{-1}
#'       \left(
#'         \mathbf{X} - \boldsymbol{\mu}
#'       \right)
#'     \right]
#'   }{
#'     \sqrt{
#'       \left(
#'         2 \pi
#'       \right)^{k}
#'       | \boldsymbol{\Sigma} |
#'     }
#'   }
#' }
#' with
#' \eqn{k}-dimensional
#' random vector
#' \eqn{\mathbf{X} \in \boldsymbol{\mu} + \textrm{span}\left(\boldsymbol{\Sigma}\right) \subseteq \mathbf{R}^k},
#' \eqn{\boldsymbol{\mu}}
#' is the location parameter mean
#' (\eqn{\boldsymbol{\mu} \in \mathbf{R}^k}),
#' and
#' \eqn{\boldsymbol{\Sigma}}
#' is the variance-covariance matrix
#' (\eqn{\boldsymbol{\Sigma} \in \mathbf{R}^{k \times k}}).
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param X Numeric matrix.
#'   Values of the \eqn{k}-dimensional random variable \eqn{X}.
#' @param mu Numeric vector.
#'   Location parameter mean vector of length \eqn{k}
#'   (\eqn{\boldsymbol{\mu}}).
#' @param Sigma Numeric matrix.
#'   \eqn{k \times k} variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma}}).
#' @importFrom stats mahalanobis
#' @references
#'   [Wikipedia: Multivariate Normal Distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
#' @export
mvnpdf <- function(X,
                   mu,
                   Sigma) {
  # Mahalanobis distance squared
  D2 <- mahalanobis(
    x = X,
    center = mu,
    cov = Sigma
  )
  (exp((-1 / 2) * D2)) / sqrt((2 * pi)^(ncol(X)) * det(Sigma))
}

#' Multivariate Normal - Log-Likelihood
#'
#' Calculates the log-likelihood of \eqn{\mathbf{X}}
#' following a mutivariate normal distribution.
#'
#' The natural log of the likelihood function for the multivariate normal
#' (or multivariate Gaussian,
#' or joint normal)
#' distribution is given by
#' \deqn{
#'   \ln \mathcal{L}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'     \mid
#'     \mathbf{X}
#'   \right)
#'   =
#'   \mathcal{l}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'     \mid
#'     \mathbf{X}
#'   \right) \\
#'   \mathcal{l}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'     \mid
#'     \mathbf{X}
#'   \right)
#'   =
#'   -
#'   \frac{1}{2}
#'   \left[
#'     \ln
#'     \left(
#'       | \boldsymbol{\Sigma} |
#'     \right)
#'     +
#'     \left(
#'       \mathbf{x} - \boldsymbol{\mu}
#'     \right)^{\prime}
#'     \boldsymbol{\Sigma}^{-1}
#'     \left(
#'       \mathbf{x} - \boldsymbol{\mu}
#'     \right)
#'     +
#'     k \ln
#'     \left(
#'     2 \pi
#'     \right)
#'   \right]
#' }
#' with
#' \eqn{k}-dimensional
#' random vector
#' \eqn{\mathbf{X} \in \boldsymbol{\mu} + \textrm{span}\left(\boldsymbol{\Sigma}\right) \subseteq \mathbf{R}^k},
#' \eqn{\boldsymbol{\mu}}
#' is the location parameter mean
#' (\eqn{\boldsymbol{\mu} \in \mathbf{R}^k}),
#' and
#' \eqn{\boldsymbol{\Sigma}}
#' is the variance-covariance matrix
#' (\eqn{\boldsymbol{\Sigma} \in \mathbf{R}^{k \times k}}).
#'
#' The negative log-likelihood is given by
#' \deqn{
#'   -
#'   \mathcal{l}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'     \mid
#'     \mathbf{X}
#'   \right)
#'   =
#'   \frac{1}{2}
#'   \left[
#'     \ln
#'     \left(
#'       | \boldsymbol{\Sigma} |
#'     \right)
#'     +
#'     \left(
#'       \mathbf{x} - \boldsymbol{\mu}
#'     \right)^{\prime}
#'     \boldsymbol{\Sigma}^{-1}
#'     \left(
#'       \mathbf{x} - \boldsymbol{\mu}
#'     \right)
#'     +
#'     k \ln
#'     \left(
#'     2 \pi
#'     \right)
#'   \right] .
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams mvnpdf
#' @inheritParams normll
#' @inherit mvnpdf references
#' @export
mvnll <- function(X,
                  mu,
                  Sigma,
                  neg = TRUE) {
  # Mahalanobis distance squared
  D2 <- mahalanobis(
    x = X,
    center = mu,
    cov = Sigma
  )
  if (neg) {
    # negative log-likelihood
    return(
      (0.5) * (log(det(Sigma)) + D2 + ncol(X) * log(2 * pi))
    )
  } else {
    # log-likelihood
    return(
      (-0.5) * (log(det(Sigma)) + D2 + ncol(X) * log(2 * pi))
    )
  }
}

#' Multivariate Normal - Two Log-Likelihood
#'
#' Calculates the two log-likelihood of \eqn{\mathbf{X}}
#' following a mutivariate normal distribution.
#'
#' The two log-likelihood for the multivariate normal
#' (or multivariate Gaussian,
#' or joint normal)
#' distribution is given by
#' \deqn{
#'   2
#'   \ln \mathcal{L}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'     \mid
#'     \mathbf{X}
#'   \right)
#'   =
#'   2
#'   \mathcal{l}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'     \mid
#'     \mathbf{X}
#'   \right) \\
#'   2
#'   \mathcal{l}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'     \mid
#'     \mathbf{X}
#'   \right)
#'   =
#'   -
#'   \ln
#'   \left(
#'     | \boldsymbol{\Sigma} |
#'   \right)
#'   -
#'   \left(
#'     \mathbf{x} - \boldsymbol{\mu}
#'   \right)^{\prime}
#'   \boldsymbol{\Sigma}^{-1}
#'   \left(
#'     \mathbf{x} - \boldsymbol{\mu}
#'   \right)
#'   -
#'   k \ln
#'   \left(
#'   2 \pi
#'   \right)
#' }
#' with
#' \eqn{k}-dimensional
#' random vector
#' \eqn{\mathbf{X} \in \boldsymbol{\mu} + \textrm{span}\left(\boldsymbol{\Sigma}\right) \subseteq \mathbf{R}^k},
#' \eqn{\boldsymbol{\mu}}
#' is the location parameter mean
#' (\eqn{\boldsymbol{\mu} \in \mathbf{R}^k}),
#' and
#' \eqn{\boldsymbol{\Sigma}}
#' is the variance-covariance matrix
#' (\eqn{\boldsymbol{\Sigma} \in \mathbf{R}^{k \times k}}).
#'
#' The negative two log-likelihood is given by
#' \deqn{
#'   -2
#'   \mathcal{l}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'     \mid
#'     \mathbf{X}
#'   \right)
#'   =
#'   \ln
#'   \left(
#'     | \boldsymbol{\Sigma} |
#'   \right)
#'   +
#'   \left(
#'     \mathbf{x} - \boldsymbol{\mu}
#'   \right)^{\prime}
#'   \boldsymbol{\Sigma}^{-1}
#'   \left(
#'     \mathbf{x} - \boldsymbol{\mu}
#'   \right)
#'   +
#'   k \ln
#'   \left(
#'   2 \pi
#'   \right) .
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams mvnpdf
#' @inheritParams norm2ll
#' @inherit mvnpdf references
#' @export
mvn2ll <- function(X,
                   mu,
                   Sigma,
                   neg = TRUE) {
  # Mahalanobis distance squared
  D2 <- mahalanobis(
    x = X,
    center = mu,
    cov = Sigma
  )
  if (neg) {
    # negative log-likelihood
    return(
      log(det(Sigma)) + D2 + ncol(X) * log(2 * pi)
    )
  } else {
    # log-likelihood
    return(
      -log(det(Sigma)) - D2 - ncol(X) * log(2 * pi)
    )
  }
}
