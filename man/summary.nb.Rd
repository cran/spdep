% Copyright 2001-3 by Roger S. Bivand
\name{summary.nb}
\alias{summary.nb}
\alias{print.nb}
\alias{summary.listw}
\alias{print.listw}
\title{Print and summary function for neighbours and weights lists}
\description{
The function prints summary measures for links in a neighbours list. If a
matrix of coordinates is given as well, summary descriptive measures for
the link lengths are also printed. Print and summary functions are also available for \code{"listw"} weights list objects.
}
\usage{
summary.nb(object, coords=NULL, lonlat = FALSE, scale = 1, ...)
print.nb(x, ...)
summary.listw(object, coords, lonlat, zero.policy = FALSE, scale = 1, ...)
print.listw(x, zero.policy = FALSE, ...)
}

\arguments{
  \item{object}{an object of class \code{nb}}
  \item{coords}{matrix of region point coordinates}
  \item{lonlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers}
  \item{...}{additional arguments affecting the output produced}
  \item{x}{an object of class \code{nb}}
  \item{zero.policy}{If FALSE stop with error for any empty neighbour sets}
  \item{scale}{passed through to \code{stem()} for control of plot length}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{plot.nb}}}

\examples{
data(columbus)
col.gal.nb
summary(col.gal.nb, coords)
col.listw <- nb2listw(col.gal.nb, style="W")
col.listw
summary(col.listw)
}
\keyword{spatial}
