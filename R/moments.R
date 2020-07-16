#' Central Moments
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
#' @author Ivan Jacob Agaloos Pesigan
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

#' Cumulant
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
#' @author Ivan Jacob Agaloos Pesigan
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

#' Kurtosis
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
#' @author Ivan Jacob Agaloos Pesigan
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
    return(out)
  } else {
    return(out + 3)
  }
}


#' Skewness
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
#' @author Ivan Jacob Agaloos Pesigan
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
