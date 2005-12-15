% Copyright 2001-4 by Roger S. Bivand
\name{poly2nb}
\alias{poly2nb}
\title{Construct neighbours list from polygon list}
\description{
The function builds a neighbours list based on regions with contiguous boundaries, that is sharing one or more boundary point. The current function is in part interpreted and may run slowly for many regions or detailed boundaries, but from 0.2-16 should not fail because of lack of memory when single polygons are built of very many border coordinates.
}
\usage{
poly2nb(pl, row.names = NULL, snap=sqrt(.Machine$double.eps), queen=TRUE)
}

\arguments{
  \item{pl}{list of polygons of class \code{polylist}, or an object of class \code{SpatialPolygons}}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}; if \code{polys} has a \code{region.id} attribute, it is copied to the neighbours list.}
  \item{snap}{boundary points less than \code{snap} distance apart are considered to indicate contiguity}
  \item{queen}{if TRUE, a single shared boundary point meets the contiguity condition, if FALSE, more than one shared point is required}
}
\value{
A neighbours list with class \code{nb}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}, \code{\link[maptools]{plot.polylist}}}

\examples{
data(columbus)
xx <- poly2nb(polys)
dxx <- diffnb(xx, col.gal.nb)
library(maptools)
plot(polys, border="grey", forcefill=FALSE)
plot(col.gal.nb, coords, add=TRUE)
plot(dxx, coords, add=TRUE, col="red")
title(main="Differences (red) in Columbus GAL weights (black)\nand polygon generated queen weights")
xxx <- poly2nb(polys, queen=FALSE)
dxxx <- diffnb(xxx, col.gal.nb)
plot(polys, border = "grey", forcefill=FALSE)
plot(col.gal.nb, coords, add = TRUE)
plot(dxxx, coords, add = TRUE, col = "red")
title(main="Differences (red) in Columbus GAL weights (black)\nand polygon generated rook weights")
cards <- card(xx)
maxconts <- which(cards == max(cards))
if(length(maxconts) > 1) maxconts <- maxconts[1]
fg <- rep("grey", length(polys))
fg[maxconts] <- "red"
fg[xx[[maxconts]]] <- "green"
plot(polys, col=fg, forcefill=FALSE)
title(main="Region with largest number of contiguities")
xx <- readShapePoly(system.file("shapes/sids.shp", package="maptools")[1], IDvar="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
xxnb <- poly2nb(as(xx, "SpatialPolygons"))
plot(xx)
plot(xxnb, coordinates(xx), add=TRUE, col="blue")
}
\keyword{spatial}
