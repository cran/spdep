% Copyright 2001 by Roger S. Bivand
\name{plotpolys}
\alias{plotpolys}
\alias{get.polylist}
\alias{convert.pl}
\title{Plot polygons}
\description{
A helper function for plotting polygons in a global bounding box. get.polylist() is a simple function to make imported GIS vector polygons into a "polylist" object from a "Map" object; convert.pl() serves to convert the deprecated "multipart" "polylist" format to an NA-separated format now used by plotpolys and poly2nb
}
\usage{
plotpolys(pl, bb, col = NA, border = par("fg"), add = FALSE, xlim, ylim)
get.polylist(Map, region.id=NULL)
convert.pl(pl)
}

\arguments{
  \item{pl}{list of polygons of class \code{polylist}}
  \item{bb}{matrix of polygon bounding boxes - columns are LL(x), LL(y), UR(x), UR(y)}
  \item{col}{colurs to use for filling the polygons}
  \item{border}{the color to draw the border}
  \item{add}{add to existing plot}
  \item{xlim, ylim}{numeric of length 2, giving the x and y coordinates ranges}
  \item{Map}{object of class \code{Map}, imported from a GIS polygon vector file using unpublished package "maptools"}
  \item{region.id}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{poly2nb}}}

\examples{
example(poly2nb)
}
\keyword{spatial}
