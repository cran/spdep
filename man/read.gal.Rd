% Copyright 2001 by Roger S. Bivand
\name{read.gal}
\alias{read.gal}
\title{Read a GAL lattice file into a neighbours list}
\description{
The function reads a GAL lattice file into a neighbours list for spatial analysis.
}
\usage{
read.gal(file, row.names=NULL)
}

\arguments{
  \item{file}{name of file with GAL lattice data}
  \item{row.names}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}, default \code{seq(1, nrow(x))}}
}
\details{
Luc Anselin (1995): SpaceStat; GAL - Geographical Algorithms
Library, University of Newcastle
}
\value{
The function returns an object of class \code{nb} with a list of integer vectors containing neighbour region number ids.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}}

\examples{
data(columbus)
data(columbus)
GALfile <- tempfile("GAL")
write.nb.gal(col.gal.nb, GALfile)
col.queen <- read.gal(GALfile)
unlink(GALfile)
summary(diffnb(col.queen, col.gal.nb))
}
\keyword{spatial}
