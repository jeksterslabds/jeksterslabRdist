#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Central Moments
#'
#' @description The sample central moment is defined by
#' \deqn{
#'   m_j
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i
#'     -
#'     \bar{x}
#'   \right)^j
#'   %(\#eq:dist-moments-sample-central)
#' }
#' where
#' - \eqn{n} is the sample size,
#' - \eqn{x = \{x_1 \dots x_n\}} is a set of values of a random variable \eqn{X},
#' - \eqn{\bar{x}} is the sample mean of \eqn{x}, and
#' - \eqn{j} is the \eqn{j}th central moment.
#'
#' @details
#' \deqn{
#'   m_0
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i
#'     -
#'     \bar{x}
#'   \right)^0
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     1
#'   \right)
#'   =
#'   \frac{1}{n}
#'   n
#'   =
#'   \frac{n}{n}
#'   =
#'   1
#'   %(\#eq:dist-moments-sample-central-zero)
#' }
#'
#' \deqn{
#'   m_1
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i
#'     -
#'     \bar{x}
#'   \right)^1
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     0
#'   \right)
#'   =
#'   \frac{1}{n}
#'   0
#'   =
#'   \frac{0}{n}
#'   =
#'   0
#'   %(\#eq:dist-moments-sample-central-first)
#' }
#'
#' \deqn{
#'   m_2
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i
#'     -
#'     \bar{x}
#'   \right)^2
#'   %(\#eq:dist-moments-sample-central-second)
#' }
#'
#' \deqn{
#'   m_3
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i
#'     -
#'     \bar{x}
#'   \right)^3
#'   %(\#eq:dist-moments-sample-central-third)
#' }
#'
#' \deqn{
#'   m_4
#'   =
#'   \frac{1}{n}
#'   \sum_{i = 1}^{n}
#'   \left(
#'     x_i
#'     -
#'     \bar{x}
#'   \right)^4
#'   %(\#eq:dist-moments-sample-central-fourth)
#' }
#'
#' - The "zeroth" central moment is 1.
#' - The first central moment is 0.
#' - The second cental moment is the variance.
#' - The third central moment is used to define skewness.
#' - The fourth central moment is used to define kurtosis.
#'
#' @family moments functions
#' @keywords moments
#' @param x Numeric vector.
#'   Sample data.
#' @param j Integer.
#'   `j`th moment.
#'   From 0 to 4.
#' @references
#'   [Wikipedia: Central Moment](https://en.wikipedia.org/wiki/Central_moment)
#'
#'   [Wikipedia: Standardized Moment](https://en.wikipedia.org/wiki/Standardized_moment)
#' @export
moment <- function(x, j) {
  if (j %in% c(0, 1, 2, 3, 4)) {
    if (j > 0 & j <= 4) {
      return((1 / length(x)) * (sum((x - mean(x))^j)))
    } else {
      stop(
        "Length of `j` should be 1."
      )
    }
  } else {
    stop(
      "Allowed values for `j` are 0, 1, 2, 3, and 4."
    )
  }
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Cumulant
#'
#' @details
#' \deqn{
#'   K_2
#'   =
#'   \frac{n}{n - 1}
#'   m_2
#' }
#'
#' \deqn{
#'   K_3
#'   =
#'   \frac{n^2}{ \left( n - 1 \right) \left( n - 2 \right) }
#'   m_3
#' }
#'
#' \deqn{
#'   K_4
#'   =
#'   \frac{n^2}{ \left( n - 1 \right) \left( n - 2 \right) \left( n - 3 \right) }
#'   \left\{
#'     \left(
#'       n + 1
#'     \right)
#'     m_4 - 3
#'     \left(
#'       n - 1
#'     \right)
#'     m_{2}^{2}
#'   \right\}
#' }
#'
#' @family moments functions
#' @keywords moments
#' @inheritParams moment
#' @export
cumulant <- function(x, j) {
  if (j %in% c(2, 3, 4)) {
    n <- length(x)
    m <- moment(x = x, j)
    if (j == 2) {
      return((n / (n - 1)) * m)
    }
    if (j == 3) {
      return((n^2 / ((n - 1) * (n - 2))) * m)
    }
    if (j == 4) {
      m2 <- moment(x = x, j = 2)
      return((n^2 / ((n - 1) * (n - 2) * (n - 3))) * ((n + 1) * m - 3 * (n - 1) * m2^2))
    }
  } else {
    stop(
      "Allowed values for `j` are 2, 3, and 4."
    )
  }
}

g1 <- function(x) {
  m2 <- moment(x = x, j = 2)
  m3 <- moment(x = x, j = 3)
  m3 / m2^(3 / 2)
}
g2 <- function(x) {
  m2 <- moment(x = x, j = 2)
  m4 <- moment(x = x, j = 4)
  (m4 / m2^2) - 3
}
G1 <- function(x) {
  K2 <- cumulant(x = x, j = 2)
  K3 <- cumulant(x = x, j = 3)
  n <- length(x)
  K3 / K2^(3 / 2)
}
G2 <- function(x) {
  K2 <- cumulant(x = x, j = 2)
  K4 <- cumulant(x = x, j = 4)
  n <- length(x)
  K4 / K2^2
}
b1 <- function(x) {
  s <- sd(x)
  m3 <- moment(x = x, j = 3)
  m3 / s^3
}
b2 <- function(x) {
  s <- sd(x)
  m4 <- moment(x = x, j = 4)
  (m4 / s^4) - 3
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Kurtosis
#'
#' @description Calculates kurtosis
#' which is a measure of the "taildness"
#' of a distribution.
#'
#' @details
#' Type 1
#' \deqn{
#'   g_2
#'   =
#'   \frac{m_4}{m_{2}^{2}}
#'   %(\#eq:dist-moments-g2)
#' }
#'
#' Type 2
#' \deqn{
#'   G_2
#'   =
#'   \frac{K_4}{K_{2}^{2}}
#'   %(\#eq:dist-moments-G2)
#' }
#'
#' Type 3
#' \deqn{
#'   b_2
#'   =
#'   \frac{m_4}{s^4}
#'   %(\#eq:dist-moments-b2)
#' }
#'
#' @family moments functions
#' @keywords moments
#' @inheritParams moment
#' @param type Integer.
#'   1, 2, or 3.
#'   See details.
#' @param excess Logical.
#'   Return excess kurtosis
#'   (kurtosis minus 3).
#' @examples
#' x <- c(
#'   10, 11, 12, 12, 12, 13, 13, 13, 13, 13,
#'   13, 14, 14, 14, 14, 15, 15, 15, 16, 17
#' )
#' kurt(x, type = 1)
#' kurt(x, type = 2)
#' kurt(x, type = 3)
#' x <- rnorm(n = 1000)
#' kurt(x, type = 1)
#' kurt(x, type = 2)
#' kurt(x, type = 3)
#' @references
#'   Joanes, D., & Gill, C. (1998).
#'   Comparing Measures of Sample Skewness and Kurtosis.
#'   *Journal of the Royal Statistical Society. Series D (The Statistician)*, *47*(1), 183-189.
#'   [2988433](https://www.jstor.org/stable/2988433)
#'
#'   [Wikipedia: Kurtosis](https://en.wikipedia.org/wiki/Kurtosis)
#' @importFrom stats sd
#' @export
kurt <- function(x, type = 2, excess = TRUE) {
  if (type == 1) {
    out <- g2(x)
  }
  if (type == 2) {
    out <- G2(x)
  }
  if (type == 3) {
    out <- b2(x)
  }
  if (excess) {
    # g2, G2, and b2 already computes excess kurtosis
    return(out)
  } else {
    return(out + 3)
  }
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Skewness
#'
#' @description Calculates kurtosis
#' which is a measure of the asymmetry
#' of a distribution.
#'
#' @details
#' Type 1
#' \deqn{
#'   g_1
#'   =
#'   \frac{m_3}{m_{2}^{\frac{3}{2}}}
#'   %(\#eq:dist-moments-g1)
#' }
#'
#' Type 2
#' \deqn{
#'   G_1
#'   =
#'   \frac{K_3}{K_{2}^{\frac{3}{2}}}
#'   %(\#eq:dist-moments-G1)
#' }
#'
#' Type 3
#' \deqn{
#'   b_1
#'   =
#'   \frac{m_3}{s^3}
#'   %(\#eq:dist-moments-b1)
#' }
#'
#' @family moments functions
#' @keywords moments
#' @inheritParams moment
#' @param type Integer.
#'   1, 2, or 3.
#'   See details.
#' @examples
#' x <- c(
#'   10, 11, 12, 12, 12, 13, 13, 13, 13, 13,
#'   13, 14, 14, 14, 14, 15, 15, 15, 16, 17
#' )
#' skew(x, type = 1)
#' skew(x, type = 2)
#' skew(x, type = 3)
#' x <- rnorm(n = 1000)
#' skew(x, type = 1)
#' skew(x, type = 2)
#' skew(x, type = 3)
#' @references
#'   Joanes, D., & Gill, C. (1998).
#'   Comparing Measures of Sample Skewness and Kurtosis.
#'   *Journal of the Royal Statistical Society. Series D (The Statistician)*, *47*(1), 183-189.
#'   [2988433](https://www.jstor.org/stable/2988433)
#'
#'   [Wikipedia: Kurtosis](https://en.wikipedia.org/wiki/Skewness)
#' @importFrom stats sd
#' @export
skew <- function(x, type = 2) {
  if (type == 1) {
    return(g1(x))
  }
  if (type == 2) {
    return(G1(x))
  }
  if (type == 3) {
    return(b1(x))
  }
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Compute mean, variance, standard deviation, skewness, and kurtosis
#'
#' @importFrom stats var sd
#' @family moments functions
#' @keywords moments
#' @inheritParams moment
#' @param type Integer.
#'   1, 2, or 3.
#'   See [`skew()`] and [`kurt()`].
#' @export
moments <- function(x, type = 2) {
  c(
    mean = mean(x = x),
    var = var(x = x),
    sd = sd(x = x),
    skew = skew(x = x, type = type),
    kurt = kurt(x = x, type = type)
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Multivariate Skewness
#'
#' @description Mardia's estimate of multivariate skewness is given by
#'   \deqn{
#'     b_{1, k}
#'     =
#'     \frac{1}{n^2}
#'     \sum_{i = 1}^{n}
#'     \sum_{j = 1}^{n}
#'     \left[
#'       \left(
#'         \mathbf{X}_{i} - \mathbf{\bar{X}}
#'       \right)^{T}
#'       \boldsymbol{\hat{\Sigma}}^{-1}
#'       \left(
#'         \mathbf{X}_{j} - \mathbf{\bar{X}}
#'       \right)
#'     \right]^{3}
#'   }
#'   where
#'   - \eqn{\mathbf{X}} is the \eqn{n \times k} sample data
#'   - \eqn{\mathbf{\bar{X}}} represent sample means
#'   - \eqn{\mathbf{X}_{i} - \mathbf{\bar{X}}} and \eqn{\mathbf{X}_{j} - \mathbf{\bar{X}}} represent deviations from the mean
#'   - \eqn{\boldsymbol{\hat{\Sigma}}} is the estimated variance-covariance matrix of \eqn{\mathbf{X}} using sample data
#'
#' @details If the null hypothesis that \eqn{\mathbf{X}} comes from a multivariate normal distribution is true,
#'   \eqn{\frac{nb_{1, k}}{6}} follows a chi-square \eqn{\left( \chi^2 \right)} distribution with a df of \eqn{\frac{k(k + 1)(k + 2)}{6}} .
#'
#' @importFrom stats pchisq cov
#' @importFrom jeksterslabRmatrix is.positive.definite
#' @param X Matrix or data frame.
#' @references
#' Mardia, K. V. (1970).
#' Measures of multivariate skewness and kurtosis with applications.
#' *Biometrika*, *57*(3), 519-530.
#' doi:[10.2307/2334770](https://doi.org/10.2307/2334770).
#'
#' Mardia, K. V. (1974).
#' Applications of Some Measures of Multivariate Skewness and Kurtosis in Testing Normality and Robustness Studies.
#' *Sankhyā: The Indian Journal of Statistics, Series B (1960-2002)*, *36*(2), 115-128.
#' @return Returns a vector with the following elements
#' \describe{
#'   \item{b1}{Estimate of multivariate skewness \eqn{\left( b_{1, k} \right)} .}
#'   \item{chisq}{chi-square statistic \eqn{\left( \frac{nb_{1, k}}{6} \right)} .}
#'   \item{df}{Degrees of freedom \eqn{\left( \frac{k(k + 1)(k + 2)}{6} \right)} .}
#'   \item{p}{p-value associated with the chi-square statistic.}
#' }
#' @examples
#' set.seed(42)
#' n <- 100
#' mu <- c(0, 0, 0)
#' Sigma <- matrix(
#'   data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1),
#'   ncol = 3
#' )
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' mardiaskew(X)
#' @export
mardiaskew <- function(X) {
  Sigma <- cov(X)
  if (!is.positive.definite(X = Sigma, stop = FALSE)) {
    stop(
      "cov(X) is non-positive definite."
    )
  }
  n <- nrow(X)
  k <- ncol(X)
  deviation <- scale(
    X,
    center = TRUE,
    scale = FALSE
  )
  invSigma <- solve(Sigma)
  # for loop version----------------------------------------------------------------
  #  deviation <- scale(
  #  X,
  #  center = TRUE,
  #  scale = FALSE
  # )
  # term <- matrix(
  #  data = NA,
  #  ncol = n,
  #  nrow = n
  # )
  # for (j in 1:n) {
  #  for (i in 1:n) {
  #    term[i, j] <- (deviation[i, ] %*% invSigma %*% as.matrix(deviation[j, ]))^3
  #  }
  # }
  # b1 <- (1 / n^2) * sum(term)
  # vectorized version-----------------------------------------------------------------
  term <- deviation %*% invSigma %*% t(deviation)
  b1 <- (1 / n^2) * sum(term^3)
  chisq <- n * (b1 / 6)
  df <- k * (k + 1) * (k + 2) / 6
  p <- pchisq(
    q = chisq,
    df = df,
    lower.tail = FALSE
  )
  c(
    b1 = b1,
    chisq = chisq,
    df = df,
    p = p
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Multivariate Kurtosis
#'
#' @description Mardia's estimate of multivariate kurtosis is given by
#'   \deqn{
#'     b_{2, k}
#'     =
#'     \frac{1}{n}
#'     \sum_{i = 1}^{n}
#'     \left[
#'       \left(
#'         \mathbf{X}_{i} - \mathbf{\bar{X}}
#'       \right)^{T}
#'       \boldsymbol{\hat{\Sigma}}^{-1}
#'       \left(
#'         \mathbf{X}_{i} - \mathbf{\bar{X}}
#'       \right)
#'     \right]^{2}
#'   }
#'   where
#'   - \eqn{\mathbf{X}} is the \eqn{n \times k} sample data
#'   - \eqn{\mathbf{\bar{X}}} represent sample means
#'   - \eqn{\mathbf{X}_{i} - \mathbf{\bar{X}}} represent deviations from the mean
#'   - \eqn{\boldsymbol{\hat{\Sigma}}} is the estimated variance-covariance matrix of \eqn{\mathbf{X}} using sample data
#'
#' @details If the null hypothesis that \eqn{\mathbf{X}} comes from a multivariate normal distribution is true,
#'   \eqn{b_{2, k}} follows a normal distribution with a mean of \eqn{k \left( k + 2 \right)}
#'   and a variance of \eqn{\frac{8 k \left( k + 2 \right)}{n}}.
#'   Consequently,
#'   \deqn{
#'     \frac{b_{2, k} - k \left( k + 2 \right)}{\sqrt{\frac{8 k \left( k + 2 \right)}{n}}}
#'   }
#'   asymptotically follows a standard normal distribution \eqn{\mathcal{N} \left( 0, 1 \right)} .
#'
#' @importFrom stats pnorm
#' @inheritParams mardiaskew
#' @inherit mardiaskew references
#' @return Returns a vector with the following elements
#' \describe{
#'   \item{b2}{Estimate of multivariate kurtosis \eqn{\left( b_{2, k} \right)} .}
#'   \item{z}{z-statistic \eqn{\left( \frac{b_{2, k} - k \left( k + 2 \right)}{\sqrt{\frac{8 k \left( k + 2 \right)}{n}}} \right)} .}
#'   \item{p}{p-value associated with the z-statistic.}
#' }
#' @examples
#' set.seed(42)
#' n <- 100
#' mu <- c(0, 0, 0)
#' Sigma <- matrix(
#'   data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1),
#'   ncol = 3
#' )
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' mardiakurt(X)
#' @export
mardiakurt <- function(X) {
  Sigma <- cov(X)
  if (!is.positive.definite(X = Sigma, stop = FALSE)) {
    stop(
      "cov(X) is non-positive definite."
    )
  }
  n <- nrow(X)
  k <- ncol(X)
  deviation <- scale(
    X,
    center = TRUE,
    scale = FALSE
  )
  invSigma <- solve(Sigma)
  # for loop version----------------------------------------------------------------
  #  deviation <- scale(
  #  X,
  #  center = TRUE,
  #  scale = FALSE
  # )
  # term <- rep(
  #  x = NA,
  #  length = n
  # )
  # for (i in 1:n) {
  #  term[i] <- (deviation[i, ] %*% invSigma %*% as.matrix(deviation[i, ]))^2
  # }
  # b2 <- (1 / n) * sum(term)
  # vectorized version-----------------------------------------------------------------
  term <- deviation %*% invSigma %*% t(deviation)
  b2 <- (1 / n) * sum(diag(term^2))
  z <- (b2 - ((k * (k + 2)))) / sqrt(8 * k * ((k + 2) / n))
  p <- 2 * pnorm(q = -abs(z))
  c(
    b2 = b2,
    z = z,
    p = p
  )
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Multivariate Skewness and Kurtosis
#' @inheritParams mardiaskew
#' @return Returns a vector with the following elements
#' \describe{
#'   \item{b1}{Estimate of multivariate skewness \eqn{\left( b_{1, k} \right)} .}
#'   \item{b1.chisq}{chi-square statistic \eqn{\left( \frac{nb_{1, k}}{6} \right)} .}
#'   \item{b1.df}{Degrees of freedom \eqn{\left( \frac{k(k + 1)(k + 2)}{6} \right)} .}
#'   \item{b1.p}{p-value associated with the chi-square statistic.}
#'   \item{b2}{Estimate of multivariate kurtosis \eqn{\left( b_{2, k} \right)} .}
#'   \item{b2.z}{z-statistic \eqn{\left( \frac{b_{2, k} - k \left( k + 2 \right)}{\sqrt{\frac{8 k \left( k + 2 \right)}{n}}} \right)} .}
#'   \item{b2.p}{p-value associated with the z-statistic.}
#' }
#' @examples
#' set.seed(42)
#' n <- 100
#' mu <- c(0, 0, 0)
#' Sigma <- matrix(
#'   data = c(1, 0.5, 0.5, 0.5, 1, 0.5, 0.5, 0.5, 1),
#'   ncol = 3
#' )
#' X <- MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
#' mardia(X)
#' @export
mardia <- function(X) {
  Sigma <- cov(X)
  if (!is.positive.definite(X = Sigma, stop = FALSE)) {
    stop(
      "cov(X) is non-positive definite."
    )
  }
  n <- nrow(X)
  k <- ncol(X)
  deviation <- scale(
    X,
    center = TRUE,
    scale = FALSE
  )
  invSigma <- solve(Sigma)
  term <- deviation %*% invSigma %*% t(deviation)
  b1 <- (1 / n^2) * sum(term^3)
  b1.chisq <- n * (b1 / 6)
  b1.df <- k * (k + 1) * (k + 2) / 6
  b1.p <- pchisq(
    q = b1.chisq,
    df = b1.df,
    lower.tail = FALSE
  )
  b2 <- (1 / n) * sum(diag(term^2))
  b2.z <- (b2 - ((k * (k + 2)))) / sqrt(8 * k * ((k + 2) / n))
  p.b2 <- 2 * pnorm(q = -abs(b2.z))
  c(
    b1 = b1,
    b1.chisq = b1.chisq,
    b1.df = b1.df,
    b1.p = b1.p,
    b2 = b2,
    b2.z = b2.z,
    p.b2 = p.b2
  )
}
