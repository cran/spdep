% Copyright 2001-3 by Roger S. Bivand
\name{poly2nb}
\alias{poly2nb}
\title{Construct neighbours list from polygon list}
\description{
The function builds a neighbours list based on regions with contiguous boundaries, that is sharing one or more boundary point. The current function is interpreted and may run slowly for many regions or detailed boundaries.
}
\usage{
poly2nb(pl, bb, row.names = NULL, snap=sqrt(.Machine$double.eps), queen=TRUE)
}

\arguments{
  \item{pl}{list of polygons of class \code{polylist}}
  \item{bb}{matrix of polygon bounding boxes - columns are LL(x), LL(y), UR(x), UR(y)}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}; if \code{polys} has a \code{region.id} attribute, it is copied to the neighbours list.}
  \item{snap}{boundary points less than \code{snap} distance apart are considered to indicate contiguity}
  \item{queen}{if TRUE, a single shared boundary point meets the contiguity condition, if FALSE, more than one shared point is required}
}
\value{
A neighbours list with class \code{nb}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}, \code{\link{plotpolys}}}

\examples{
data(columbus)
xx <- poly2nb(polys, bbs)
dxx <- diffnb(xx, col.gal.nb)
plotpolys(polys, bbs, border="grey")
plot(col.gal.nb, coords, add=TRUE)
plot(dxx, coords, add=TRUE, col="red")
title(main="Differences (red) in Columbus GAL weights (black)\nand polygon generated queen weights")
xxx <- poly2nb(polys, bbs, queen=FALSE)
dxxx <- diffnb(xxx, col.gal.nb)
plotpolys(polys, bbs, border = "grey")
plot(col.gal.nb, coords, add = TRUE)
plot(dxxx, coords, add = TRUE, col = "red")
title(main="Differences (red) in Columbus GAL weights (black)\nand polygon generated rook weights")
cards <- card(xx)
maxconts <- which(cards == max(cards))
if(length(maxconts) > 1) maxconts <- maxconts[1]
fg <- rep("grey", length(polys))
fg[maxconts] <- "red"
fg[xx[[maxconts]]] <- "green"
plotpolys(polys, bbs, col=fg)
title(main="Region with largest number of contiguities")
}
\keyword{spatial}
