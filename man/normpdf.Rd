% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{normpdf}
\alias{normpdf}
\title{Normal Distribution - Probablity Density Function}
\usage{
normpdf(x, mu = 0, sigma = 1, log = FALSE)
}
\arguments{
\item{x}{Numeric vector.
Values of the random variable \eqn{X}.}

\item{mu}{Numeric.
Location parameter mean \eqn{\mu}.}

\item{sigma}{Numeric.
Positive number.
Scale parameter standard deviation \eqn{\sigma = \sqrt{\sigma^2}}.}

\item{log}{Logical.
If \code{TRUE},
returns \eqn{\log \left( f \left( x \right) \right)}.}
}
\value{
Returns \eqn{f \left( x \right)}
using the probablity density function
with the supplied parameter/s.
If \code{log = TRUE},
returns \eqn{\log \left( f \left( x \right) \right)}.
}
\description{
Calculates probablities
from the probability density function
of the normal distribution
\eqn{
  X
  \sim
  \mathcal{N}
  \left(
    \mu,
    \sigma^2
  \right)
  %(\#eq:dist-X-norm)
} .
This function is identical to \code{\link[=dnorm]{dnorm()}}.
}
\details{
The normal (or Gaussian or Gauss or Laplace–Gauss)
distribution is given by
\deqn{
  X
  \sim
  \mathcal{N}
  \left(
    \mu,
    \sigma^2
  \right)
  %(\#eq:dist-X-norm)
}
and has
the probability density function (PDF)
\deqn{
  f
  \left(
    x
  \right)
  =
  \frac{1}{\sigma \sqrt{2 \pi}}
  \exp
  \left[
    -
    \frac{1}{2}
    \left(
      \frac{x - \mu}{\sigma}
    \right)^2
  \right]
  %(\#eq:dist-normpdf-1)
}
or
\deqn{
  f
  \left(
    x
  \right)
  =
  \frac{1}{\sqrt{2 \pi \sigma^2}}
  \exp
    \left[
      -
      \frac{
      \left(
        x - \mu
      \right)^2}
      {2 \sigma^2}
    \right]
  %(\#eq:dist-normpdf-2)
}
with
\itemize{
\item \eqn{x \in \mathbf{R}},
\item \eqn{\mu} is the location parameter mean
\eqn{\left( \mu \in \mathbf{R} \right)}, and
\item \eqn{\sigma^2} is the scale parameter variance
\eqn{\left( \sigma^2 > 0 \right)}.
}
}
\references{
\href{https://en.wikipedia.org/wiki/Normal_distribution}{Wikipedia: Normal Distribution}
}
\seealso{
Other normal likelihood functions: 
\code{\link{norm2ll}()},
\code{\link{normL}()},
\code{\link{normll}()},
\code{\link{normobj}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{normal likelihood functions}
\keyword{normal}
