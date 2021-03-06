\name{ipacose}
\alias{ipacose}
\title{The iPACOSE algorithm for graph estimation}
\description{The iPACOSE algorithm allowing to iteratively estimate a graph from a dataset and a first estimation of the partial correlation matrix.}
\usage{
ipacose(x = x, pc = pc, method = "pacose.ridge", cutoff, gr=NULL, k = 2, cv.method = "CV" ,adaptive=FALSE)
}
\arguments{
  \item{x}{a dataset (matrix) of dimensions n x p.}
  \item{pc}{a first version of the  p x p partial correlation matrix.}
  \item{method}{should be "pacose.pls", "pacose.ridge" or "pacose.adalasso".}
  \item{cutoff}{a threshold to apply to the new partial correlation matrices obtained with PACOSE in order to obtain a new version of the graph.}
  \item{gr}{the true graph (an object of class \code{\link[igraph:igraph-package]{igraph}}), when it is known. It must be set to NULL (the default value) when the true graph is not known.}
  \item{k}{integer, the number of folds to be used when selecting the optimal value for he regularization parameter.}
  \item{adaptive}{boolean, default to FALSE.} 
  \item{cv.method}{only valid when the reference method is \code{ridge.net}: determines the way the ridge parameter is computed ("HKB" or "CV", see \link{pacose.ridge}).  }
}
\value{
\item{Niter}{number of iterations.}
\item{gr_it}{the final estimated graph.}
\item{pcor_it}{the final estimated p x p partial correlation matrix.}
\item{sens}{when the true graph is known, the sensitivity of the estimated graph.}
\item{ppv}{when the true graph is known, the PPV of the estimated graph.}
}
\references{Guillemot V., Bender A., Boulesteix A.-L. (2012). Iterative reconstruction of high-dimensional Gaussian graphical models based on a new method to estimate partial correlations under constraints. Submitted.}
\author{Vincent Guillemot}
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
gr <- graph.adjacency((omega!=0),mode="undirected",diag=FALSE)
x <- rmvnorm(n=100,sigma=solve(omega))

# First estimation with ridge.net, the threshold is, for this toy example,
# arbitrarily set to 0.05
pcor0.0 <- ridge.net(x)$pcor
pcor0 <- (abs(pcor0.0) > 0.05)+0

# Use iPACOSE to estimate iteratively the partial correlation matrix
pcor1 <- ipacose(x=x,pc=pcor0,method="pacose.ridge",
                 cutoff=0.05,gr=NULL,cv.method="HKB")
cat("Number of iterations: ", pcor1$Niter)

invcov2pcor(omega)
pcor1$pcor_it
}
\keyword{ algebra }
\keyword{ multivariate }
