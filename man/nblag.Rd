% Copyright 2001 by Roger S. Bivand
\name{nblag}
\alias{nblag}
\title{Higher order neighbours lists}
\description{
The function creates higher order neighbour lists, where higher order neighbours are only \code{lags} links from each other on the graph described by the input neighbours list. It will refuse to lag neighbours lists with the attribute self.included set to TRUE.
}
\usage{
nblag(neighbours, maxlag)
}
\arguments{
  \item{neighbours}{input neighbours list of class \code{nb}}
  \item{maxlag}{the maximum lag to be constructed}
}
\value{
returns a list of lagged neighbours lists each with class \code{nb}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}}

\examples{
data(columbus)
summary(col.gal.nb, coords)
library(maptools)
plot(polys, border="grey", forcefill=FALSE)
plot(col.gal.nb, coords, add=TRUE)
title(main="GAL order 1 (black) and 2 (red) links")
col.lags <- nblag(col.gal.nb, 2)
summary(col.lags[[2]], coords)
plot(col.lags[[2]], coords, add=TRUE, col="red", lty=2)
}
\keyword{spatial}
