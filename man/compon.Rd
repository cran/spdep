% Copyright 2001 by Roger S. Bivand and Nicholas Lewin-Koh
\name{Graph Components}
\alias{n.comp.nb}
%\alias{reach.ij}
%\alias{strong.comp}
\title{Depth First Search on Neighbor Lists}
\description{
  \code{n.comp.nb()} finds the number of disjoint connected subgraphs in the graph depicted by \code{nb.obj} - a spatial neighbours list object.
}
\usage{
n.comp.nb(nb.obj)
}
\arguments{
 \item{nb.obj}{a neighbours list object of class \code{nb}}
}
\value{
A list of:
\item{nc}{number of disjoint connected subgraphs}
\item{comp.id}{vector with the indices of the disjoint connected subgraphs that
 the nodes in \code{nb.obj} belong to}
}
\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}
\seealso{\code{\link{plot.nb}}}
\examples{
example(columbus)
coords <- coordinates(columbus)
plot(col.gal.nb, coords, col="grey")
col2 <- droplinks(col.gal.nb, 21)
plot(col2, coords, add=TRUE)
res <- n.comp.nb(col2)
table(res$comp.id)
points(coords, col=res$comp.id, pch=16)
}

\keyword{spatial}
