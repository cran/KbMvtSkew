% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/KbSkewness.R
\name{BowleySkew}
\alias{BowleySkew}
\title{Bowley's Univariate Skewness}
\usage{
BowleySkew(x)
}
\arguments{
\item{x}{a vector of original observations.}
}
\value{
\code{BowleySkew} gives the Bowley's univariate skewness of the data.
}
\description{
Compute Bowley's Univariate Skewness.
}
\details{
Bowley's skewness is defined in terms of quantiles as
\deqn{\hat{\gamma} = \frac{Q_3 + Q_1 - 2 Q_2}{Q_3 - Q_1}}
where \eqn{Q_i} is the \eqn{i}th quartile \eqn{i=1,2,3} of the data.
}
\examples{
# Compute Bowley's univariate skewness

set.seed(2019)
x <- rnorm(1000) # Normal Distribution
BowleySkew(x)

set.seed(2019)
y <- rlnorm(1000, meanlog = 1, sdlog = 0.25) # Log-normal Distribution
BowleySkew(y)

}
\references{
Bowley, A. L. (1920). \emph{Elements of Statistics}. London : P.S. King & Son, Ltd.
}
