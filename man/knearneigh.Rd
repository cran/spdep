% Copyright 2001-3 by Roger S. Bivand
\name{knearneigh}
\alias{knearneigh}
\title{K nearest neighbours for spatial weights}
\description{
The function returns a matrix with the indices of regions belonging to the set of the k nearest neighbours of each other. If lonlat = TRUE, Great Circle distances are used.
}
\usage{
knearneigh(x, k=1, lonlat = FALSE)
}

\arguments{
  \item{x}{matrix of region point coordinates}
  \item{k}{number of nearest neighbours to be returned}
  \item{lonlat}{TRUE if point coordinates are longitude-latitude decimal degrees, in which case distances are measured in kilometers}
}
\details{
The underlying C code is based on the \code{knn} function in the class package in the VR bundle.
}
\value{
A list of class \code{knn}
  \item{nn}{integer matrix of region number ids}
  \item{np}{number of input points}
  \item{k}{input required k}
  \item{dimension}{number of columns of x}
  \item{x}{input coordinates}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{knn}}, \code{\link{dnearneigh}},
\code{\link{knn2nb}}}

\examples{
data(columbus)
col.knn <- knearneigh(coords, k=4)
plotpolys(polys, bbs, border="grey")
plot(knn2nb(col.knn), coords, add=TRUE)
title(main="K nearest neighbours, k = 4")
data(state)
us48.fipsno <- read.geoda(system.file("etc/weights/us48.txt",
 package="spdep")[1])
m50.48 <- match(us48.fipsno$State.name, state.name)
xy <- as.matrix(as.data.frame(state.center))[m50.48,]
llk4.nb <- knn2nb(knearneigh(xy, k=4, lonlat=FALSE))
gck4.nb <- knn2nb(knearneigh(xy, k=4, lonlat=TRUE))
plot(llk4.nb, xy)
plot(diffnb(llk4.nb, gck4.nb), xy, add=TRUE, col="red", lty=2)
title(main="Differences between Euclidean and Great Circle k=4 neighbours")
summary(llk4.nb, xy, lonlat=TRUE)
summary(gck4.nb, xy, lonlat=TRUE)
}
\keyword{spatial}
