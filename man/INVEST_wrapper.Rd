\name{INVEST_wrapper}
\alias{INVEST_wrapper}
\title{Wrapper for the C implementation of Wermuth's covariance selection algorithm}
\description{An R wrapper for a C implementation of Wermuth's covariance selection algorithm. Takes as an input a dataset and the indices where there should be zeros in the precision matrix and returns the corresponding precision matrix (if possible) according to Wermuth's algorithm.}
\usage{INVEST_wrapper(x, ind, delta = 1e-09, itermax = 2000)}
\arguments{
  \item{x}{a dataset (matrix) of dimensions n x p.}
  \item{ind}{matrix of the indices corresponding to the zeros in the inverse covariance matrix.}
  \item{delta}{positive real, value under which a coefficient in the inverse covariance matrix is considered to be equal to 0.}
  \item{itermax}{integer, a maximum number of iterations to avoid infinite loops.}
}
\value{The desired inverse covariance matrix with zeros where specified.}
\references{ Wermuth, N. and Scheidt, E. (1977). Fitting a covariance selection model to a matrix, algorithm 105. Journal of the Royal Statistical Society C, 26:88--92.}
\author{Andreas Bender, Vincent Guillemot}
\keyword{ algebra }
\keyword{ multivariate }