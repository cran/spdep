% Copyright 2001 by Roger S. Bivand
\name{dnearneigh}
\alias{dnearneigh}
\title{Neighbourhood contiguity by distance}
\description{
The function identifies neighbours of region points by Euclidean distance between lower (greater than) and upper (less than or equal to) bounds.
}
\usage{
dnearneigh(x, d1, d2, row.names = NULL)
}
\arguments{
  \item{x}{matrix of point coordinates}
  \item{d1}{lower distance bound}
  \item{d2}{upper distance bound}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
}
\value{
The function returns a list of integer vectors giving the region id numbers
for neighbours satisfying the distance criteria.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{knearneigh}}}

\examples{
data(columbus)
k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
col.nb.0.all <- dnearneigh(coords, 0, all.linked, row.names=rownames(columbus))
summary(col.nb.0.all, coords)
plotpolys(polys, bbs, border="grey")
plot(col.nb.0.all, coords, add=TRUE)
title(main=paste("Distance based neighbours 0-",  format(all.linked),
 " distance units", sep=""))
}
\keyword{spatial}
