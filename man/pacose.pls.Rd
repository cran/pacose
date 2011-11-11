\name{pacose.pls}
\alias{pacose.pls}
\title{The PLS version of the PACOSE algorithm}
\description{The function takes as an argument a dataset and a graph and returns an estimation of the partial correlation matrix.}
\usage{
pacose.pls(X, gg, scale = TRUE, k = 10, Ncomp = NULL, verbose = FALSE, cv.method = "CV")
}
\arguments{
  \item{X}{a dataset (matrix) of dimensions n x p.}
  \item{gg}{the graph (an object of class \code{\link[igraph0:igraph-package]{igraph0}}) to integrate.}
  \item{scale}{boolean, whether to scale the data or not, default to FALSE.}
  \item{k}{integer, the number of folds to be used when selecting the optimal value among \code{Ncomp}.}
  \item{Ncomp}{integer, the maximal number of PLS components to be used.}
  \item{verbose}{boolean, whether to print out intermediate messages or not, default to FALSE.}
  \item{cv.method}{equals "CV", for a cross validation determination of the regularization parameter. Alternative values are only used in function \code{\link{pacose.ridge}}.}
}
\value{
\item{pcor}{an estimate of the p x p partial correlation matrix.}
\item{invcov }{an estimate of the p x p precision matrix.}
\item{m }{vector containing all the optimal values for the regularization parameters.}
}
\references{Guillemot V., Bender A., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.}
\author{Vincent Guillemot}
\seealso{\code{\link[parcor:pls.net]{pls.net}}, \code{\link{pacose.ridge}}, \code{\link{pacose.adalasso}}}
\examples{
require(mvtnorm)
require(igraph0)

omega <- matrix(c(1     , -0.477, 0.304, 0.478, -0.591, 0    ,
               -0.477, 2     , 0.206, 0    , 0.382 , 0    ,
               0.304 , 0.206 , 1    , 0    , 0.181 , 0.242,
               0.478 , 0     , 0    , 3    , 0.141 , 0    ,
               -0.591, 0.382 , 0.181, 0.141, 1     , 0    ,
               0     , 0     , 0.242, 0    , 0     , 2    ),
             nrow = 6, ncol = 6)
gr <- graph.adjacency((omega!=0),mode="undirected",diag=FALSE)
x <- rmvnorm(n=100,sigma=solve(omega))
omega.hat <- pacose.pls(x,gr,verbose=TRUE)$invcov

omega
round(omega.hat,3)
}
\keyword{ algebra }
\keyword{ multivariate }
