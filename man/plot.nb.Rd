% Copyright 2001 by Roger S. Bivand
\name{plot.nb}
\alias{plot.nb}
\title{Plot a neighbours list}
\description{
A function to plot a neighbours list given point coordinates to represent the region in two dimensions.
}
\usage{
plot.nb(x, coords, col="black", points=TRUE, add=FALSE, ...)
}

\arguments{
  \item{x}{an object of class \code{nb}}
  \item{coords}{matrix of region point coordinates}
  \item{col}{plotting colour}
  \item{points}{(logical) add points to plot}
  \item{add}{(logical) add to existing plot}
  \item{...}{further graphical parameters as in \code{par(..)}}
}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}}

\examples{
data(columbus)
plot(col.gal.nb, coords)
title(main="GAL order 1 links with first nearest neighbours in red")
col.knn <- knearneigh(coords, k=1)
plot(knn2nb(col.knn), coords, add=TRUE, col="red", lty=2)
}
\keyword{spatial}
