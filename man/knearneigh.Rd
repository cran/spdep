% Copyright 2001 by Roger S. Bivand
\name{knearneigh}
\alias{knearneigh}
\title{K nearest neighbours for spatial weights}
\description{
The function returns a matrix with the indices of regions belonging to the set of the k nearest neighbours of each other.
}
\usage{
knearneigh(x, k=1)
}

\arguments{
  \item{x}{matrix of region point coordinates}
  \item{k}{number of nearest neighbours to be returned}
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
}
\keyword{spatial}
