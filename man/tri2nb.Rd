% Copyright 2001 by Roger S. Bivand
\name{tri2nb}
\alias{tri2nb}
\title{Neighbours list from tri object}
\description{
The function uses the \code{tripack} package to convert a matrix of
two-dimensional coordinates
into a neighbours list of class \code{nb} with a list of integer vectors 
containing neighbour region number ids.
}
\usage{
tri2nb(coords, row.names = NULL)
}
\arguments{
  \item{coords}{matrix of point coordinates with two columns}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
}
\value{
The function returns an object of class \code{nb} with a list of integer vectors containing neighbour region number ids.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{knn2nb}}}

\examples{
require(tripack)
data(columbus)
col.tri.nb <- tri2nb(coords, row.names=rownames(columbus))
library(maptools)
plot(polys, border="grey")
plot(col.tri.nb, coords, add=TRUE)
title(main="Raw triangulation links")
}
\keyword{spatial}
