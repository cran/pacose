\name{PACOSE-internal}
\title{Internal PACOSE functions}
\alias{Beta2parcor}
\alias{adalasso}
\alias{mylars}
\alias{nbr}
\alias{nbri}
\alias{nbru}
\alias{oneK}
\alias{oneK_MVUE}
\alias{oneK_SURE}
\alias{oppdiag}
\alias{ridge.cv}
\alias{ridge.net}
\alias{tr}
\alias{trk}
\alias{wermuthC}
\description{Internal PACOSE functions}
\usage{
Beta2parcor(Beta, verbose = FALSE)
adalasso(X, y, k = 10, use.Gram = TRUE, both = TRUE)
mylars(X, y, k = 10, use.Gram = TRUE, normalize = TRUE)
nbr(gg)
nbri(gg1, gg2)
nbru(gg1, gg2)
oneK(ind, x)
oneK_MVUE(ind, x)
oneK_SURE(ind, x, d)
oppdiag(mat)
ridge.cv(X, y, lambda = NULL, scale = TRUE, k = 10, plot.it =FALSE)
ridge.net(X, lambda = NULL, plot.it = FALSE, scale = TRUE, k =10, verbose = FALSE, cv.method = "HKB")
tr(mat)
trk(ind, x)
}
\author{Vincent Guillemot}
\keyword{internal}




