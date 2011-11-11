\name{pacose.adalasso}
\alias{pacose.adalasso}
\title{The LASSO and adaptive LASSO versions of the PACOSE algorithm}
\description{The function takes as an argument a dataset and a graph and returns an estimation of the partial correlation matrix.}
\usage{
pacose.adalasso(X, gg, k = 10, use.Gram = FALSE, both = TRUE, verbose = FALSE, cv.method = "CV")
}
\arguments{
  \item{X}{a dataset (matrix) of dimensions n x p.}
  \item{gg}{the graph (an object of class \code{\link[igraph0:igraph-package]{igraph0}}) to integrate.}
  \item{k}{integer, the number of folds to be used when selecting the optimal value of the regularization parameter.}
  \item{use.Gram}{see \code{\link[parcor:adalasso]{adalasso}}}
  \item{both}{boolean, whether to compute both the adaptive LASSO and the non adaptive LASSO versions of the partial correlation matrix or not, default to TRUE.
}
  \item{verbose}{boolean, whether to print out intermediate messages or not, default to FALSE.}
  \item{cv.method}{equals "CV", for a cross validation determination of the regularization parameter. Alternative values are only used in function \code{\link{pacose.ridge}}.}
}
\value{
list(pcor.lasso = pcor.lasso, pcor.adalasso = pcor.adalasso,
                invcov.lasso = invcov.lasso, invcov.adalasso = invcov.adalasso)
\item{pcor.lasso}{an estimate of the p x p partial correlation matrix with the non adaptive version of the method LASSO.}
\item{invcov.lasso }{an estimate of the p x p precision matrix with the non adaptive version of the method LASSO.}
\item{pcor.adalasso}{an estimate of the p x p partial correlation matrix with the adaptive version of the method LASSO.}
\item{invcov.adalasso }{an estimate of the p x p precision matrix with the adaptive version of the method LASSO.}
}
\references{Guillemot V., Bender A., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.}
\author{Vincent Guillemot}
\seealso{\code{\link[parcor:adalasso.net]{adalasso.net}}, \code{\link{pacose.ridge}}, \code{\link{pacose.pls}}}
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
omega.hat <- pacose.adalasso(x,gr)$invcov.lasso

omega
round(omega.hat,3)
}
\keyword{ algebra }
\keyword{ multivariate }
