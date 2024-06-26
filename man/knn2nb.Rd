% Copyright 2001 by Roger S. Bivand
\name{knn2nb}
\alias{knn2nb}
\title{Neighbours list from knn object}
\description{
The function converts a \code{knn} object returned by \code{knearneigh} 
into a neighbours list of class \code{nb} with a list of integer vectors 
containing neighbour region number ids.
}
\usage{
knn2nb(knn, row.names = NULL, sym = FALSE)
}

\arguments{
  \item{knn}{A knn object returned by \code{knearneigh}}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
  \item{sym}{force the output neighbours list to symmetry}
}
\value{
The function returns an object of class \code{nb} with a list of integer vectors containing neighbour region number ids. See \code{\link{card}} for details of \dQuote{nb} objects.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{knearneigh}}, \code{\link{card}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
coords <- st_coordinates(st_centroid(columbus))
col.knn <- knearneigh(coords, k=4)
plot(st_geometry(columbus), border="grey")
plot(knn2nb(col.knn), coords, add=TRUE)
title(main="K nearest neighbours, k = 4")
# example of reading points with readr::read_csv() yielding a tibble
load(system.file("etc/misc/coords.rda", package="spdep"))
class(coords)
knn2nb(knearneigh(coords, k=4))
}
\keyword{spatial}
