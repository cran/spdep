% Copyright 2001-3 by Roger S. Bivand
\name{dnearneigh}
\alias{dnearneigh}
\title{Neighbourhood contiguity by distance}
\description{
The function identifies neighbours of region points by Euclidean distance between lower (greater than) and upper (less than or equal to) bounds, or with lonlat = TRUE, by Great Circle distance in kilometers.
}
\usage{
dnearneigh(x, d1, d2, row.names = NULL, lonlat = FALSE)
}
\arguments{
  \item{x}{matrix of point coordinates}
  \item{d1}{lower distance bound}
  \item{d2}{upper distance bound}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
  \item{lonlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers}
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
library(maptools)
plot(polys, border="grey")
plot(col.nb.0.all, coords, add=TRUE)
title(main=paste("Distance based neighbours 0-",  format(all.linked),
 " distance units", sep=""))
data(state)
us48.fipsno <- read.geoda(system.file("etc/weights/us48.txt",
 package="spdep")[1])
if (as.numeric(paste(version$major, version$minor, sep="")) < 19) {
 m50.48 <- match(us48.fipsno$"State.name", state.name)
} else {
 m50.48 <- match(us48.fipsno$"State_name", state.name)
}
xy <- as.matrix(as.data.frame(state.center))[m50.48,]
llk1 <- knn2nb(knearneigh(xy, k=1, lonlat=FALSE))
all.linked <- max(unlist(nbdists(llk1, xy, lonlat=FALSE)))
ll.nb <- dnearneigh(xy, 0, all.linked, lonlat=FALSE)
summary(ll.nb, xy, lonlat=TRUE, scale=0.5)
gck1 <- knn2nb(knearneigh(xy, k=1, lonlat=TRUE))
all.linked <- max(unlist(nbdists(gck1, xy, lonlat=TRUE)))
gc.nb <- dnearneigh(xy, 0, all.linked, lonlat=TRUE)
summary(gc.nb, xy, lonlat=TRUE, scale=0.5)
plot(ll.nb, xy)
plot(diffnb(ll.nb, gc.nb), xy, add=TRUE, col="red", lty=2)
title(main="Differences between Euclidean and Great Circle neighbours")
}
\keyword{spatial}
