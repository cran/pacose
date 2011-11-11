\name{pcor2invcov}
\alias{pcor2invcov}
\title{Transforms a partial correlation matrix into an inverse covariance matrix}
\description{Transforms an input p x p positive definite partial correlation matrix into an inverse covariance matrix given the vector of the partial variances.}
\usage{pcor2invcov(A, pv)}
\arguments{
  \item{A}{a given partial correlation matrix.}
  \item{pv}{a vector of partial variances.}
}
\value{The inverse covariance matrix deduced from the given partial correlations and variances.}
\references{Whittaker, J. (1990). Graphical models in applied multivariate statistics. Wiley.}
\keyword{ algebra }
\keyword{ multivariate }

