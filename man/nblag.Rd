% Copyright 2001 by Roger S. Bivand
\name{nblag}
\alias{nblag}
\alias{nblag_cumul}
\title{Higher order neighbours lists}
\description{
The function creates higher order neighbour lists, where higher order neighbours are only \code{lags} links from each other on the graph described by the input neighbours list. It will refuse to lag neighbours lists with the attribute self.included set to TRUE. \code{nblag_cumul} cumulates neighbour lists to a single neighbour list (\dQuote{nb} object).
}
\usage{
nblag(neighbours, maxlag)
nblag_cumul(nblags)
}
\arguments{
  \item{neighbours}{input neighbours list of class \code{nb}}
  \item{maxlag}{the maximum lag to be constructed}
  \item{nblags}{a list of neighbour lists as output by code{nblag}}
}
\value{
returns a list of lagged neighbours lists each with class \code{nb}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Giovanni Millo}

\seealso{\code{\link{summary.nb}}}

\examples{
example(columbus)
coords <- coordinates(columbus)
summary(col.gal.nb, coords)
plot(columbus, border="grey")
plot(col.gal.nb, coords, add=TRUE)
title(main="GAL order 1 (black) and 2 (red) links")
col.lags <- nblag(col.gal.nb, 2)
lapply(col.lags, print)
summary(col.lags[[2]], coords)
plot(col.lags[[2]], coords, add=TRUE, col="red", lty=2)
cuml <- nblag_cumul(col.lags)
cuml
}
\keyword{spatial}
