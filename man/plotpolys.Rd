% Copyright 2001-3 by Roger S. Bivand
\name{plotpolys}
\alias{plotpolys}
\alias{leglabs}
\alias{findInterval2}
\title{Plot polygons}
\description{
A helper function for plotting polygons in a global bounding box. \code{findInterval2} is a modified version of \code{findInterval2} that permits the vector of break points to start at -Inf and/or finish at Inf. \code{leglabs} makes character strings from the same break points.
}
\usage{
plotpolys(pl, bb, col = NA, border = par("fg"), add = FALSE, xlim, ylim, ...)
leglabs(vec, under="under", over="over", between="-")
findInterval2(x, vec, rightmost.closed = FALSE, all.inside = TRUE)
}

\arguments{
  \item{pl}{list of polygons of class \code{polylist}}
  \item{bb}{matrix of polygon bounding boxes - columns are LL(x), LL(y), UR(x), UR(y)}
  \item{col}{colurs to use for filling the polygons}
  \item{border}{the color to draw the border}
  \item{add}{add to existing plot}
  \item{xlim, ylim}{numeric of length 2, giving the x and y coordinates ranges}
  \item{...}{other arguments passed to polygon}
  \item{vec}{vector of break values}
  \item{under}{character value for under}
  \item{over}{character value for over}
  \item{between}{character value for between}
  \item{x}{numeric vector to be broken into class intervals}
  \item{rightmost.closed}{as in findInterval()}
  \item{all.inside}{as in findInterval()}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{Map2poly}}, \code{\link{findInterval}}}

\examples{
data(columbus)
brks <- round(quantile(columbus$CRIME, probs=seq(0,1,0.2)), digits=2)
colours <- c("salmon1", "salmon2", "red3", "brown", "black")
plotpolys(polys, bbs, col=colours[findInterval2(columbus$CRIME, brks)])
legend(x=c(5.8, 7.1), y=c(13, 14.5), legend=leglabs(brks),
  fill=colours, bty="n")
title(main=paste("Columbus OH: residential burglaries and vehicle\nthefts",
  "per thousand households, 1980"))
}
\keyword{spatial}
