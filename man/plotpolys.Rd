% Copyright 2001 by Roger S. Bivand
\name{plotpolys}
\alias{plotpolys}
\title{Plot polygons}
\description{
A helper function for plotting polygons in a global bounding box.
}
\usage{
plotpolys(pl, bb, col = NA, border = par("fg"), add = FALSE, xlim, ylim)
}

\arguments{
  \item{pl}{list of polygons og class \code{polylist}}
  \item{bb}{matrix of polygon bounding boxes - columns are LL(x), LL(y), UR(x), UR(y)}
  \item{col}{colurs to use for filling the polygons}
  \item{border}{the color to draw the border}
  \item{add}{add to existing plot}
  \item{xlim, ylim}{numeric of length 2, giving the x and y coordinates ranges}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{poly2nb}}}

\examples{
example(poly2nb)
}
\keyword{spatial}
