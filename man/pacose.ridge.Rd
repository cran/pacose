\name{pacose.ridge}
\alias{pacose.ridge}
\title{The Ridge version of the PACOSE algorithm}
\description{The function takes as an argument a dataset and a graph and returns an estimation of the partial correlation matrix.}
\usage{
pacose.ridge(X, gg, lambda = NULL, plot.it = FALSE, scale = FALSE, k = 10, verbose = FALSE, cv.method = "CV")
}
\arguments{
  \item{X}{a dataset (matrix) of dimensions n x p.}
  \item{gg}{the graph (an object of class \code{\link[igraph:igraph-package]{igraph}}) to integrate.}
  \item{lambda}{a set (vector) of values among which to choose when estimating the partial correlations.}
  \item{plot.it}{boolean, whether to generate a graphical output or not, default to FALSE.}
  \item{scale}{boolean, whether to scale the data or not, default to FALSE.}
  \item{k}{integer, the number of folds to be used when selecting the optimal value among \code{lambda}.}
  \item{verbose}{boolean, whether to print out intermediate messages or not, default to FALSE.}
  \item{cv.method}{determines the way the ridge parameter is computed: should be equal to either "HKB" for an analytical determination, or "CV", for a cross validation alternative. When equal to "HKB", neither \code{lambda} nor \code{k} are used.}
}
\value{
\item{pcor}{an estimate of the p x p partial correlation matrix.}
\item{invcov }{an estimate of the p x p precision matrix.}
}
\references{Guillemot V., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.}
\author{Vincent Guillemot}
\seealso{\code{\link[parcor:ridge.net]{ridge.net}}, \code{\link{pacose.pls}}, \code{\link{pacose.adalasso}}}
\examples{
require(mvtnorm)
require(igraph)

omega <- matrix(c(1     , -0.477, 0.304, 0.478, -0.591, 0    ,
               -0.477, 2     , 0.206, 0    , 0.382 , 0    ,
               0.304 , 0.206 , 1    , 0    , 0.181 , 0.242,
               0.478 , 0     , 0    , 3    , 0.141 , 0    ,
               -0.591, 0.382 , 0.181, 0.141, 1     , 0    ,
               0     , 0     , 0.242, 0    , 0     , 2    ),
             nrow = 6, ncol = 6)
gr <- graph.adjacency(omega!=0,mode="undirected",diag=FALSE)
x <- rmvnorm(n=100,sigma=solve(omega))
omega.hat <- pacose.ridge(x,gr,cv.method="HKB")$invcov

omega
round(omega.hat,3)
}
\keyword{ algebra }
\keyword{ multivariate }
