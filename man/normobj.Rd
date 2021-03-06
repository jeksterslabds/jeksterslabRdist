% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/norm.R
\name{normobj}
\alias{normobj}
\title{Normal Distribution - Objective Function}
\usage{
normobj(theta, x, neg = TRUE)
}
\arguments{
\item{theta}{Vector of parameters \eqn{\theta} of the normal distribution
(\code{theta[1] = mu} (\eqn{\mu}) and \code{theta[2] = sigma} (\eqn{\sigma})).}

\item{x}{Numeric vector.
Values of the random variable \eqn{X}.}

\item{neg}{Logical.
If \code{TRUE},
returns,
negative log-likelihood.}
}
\description{
Objective function to minimize/maximize
to estimate parameters of the normal distribution.
See \code{\link[=norm2ll]{norm2ll()}}.
}
\references{
\href{https://en.wikipedia.org/wiki/normaldistribution}{Wikipedia: Normal Distribution}

\href{https://en.wikipedia.org/wiki/Independent_and_identically_distributed_random_variables}{Wikipedia: IID}

\href{https://en.wikipedia.org/wiki/Likelihood_function}{Wikipedia: Likelihood Function}
}
\seealso{
Other normal likelihood functions: 
\code{\link{norm2ll}()},
\code{\link{normL}()},
\code{\link{normll}()},
\code{\link{normpdf}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{normal likelihood functions}
\keyword{normal}
