% Copyright 2001 by Roger S. Bivand
\name{write.nb.gal}
\alias{write.nb.gal}

\title{Write a neighbours list as a GAL lattice file}
\description{
  
}
\usage{
write.nb.gal(nb, file)
}

\arguments{
  \item{nb}{an object of class \code{nb} with a list of integer vectors containing neighbour region number ids.}
  \item{file}{name of file with GAL lattice data}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{read.gal}}}

\examples{
data(columbus)
GALfile <- tempfile("GAL")
write.nb.gal(col.gal.nb, GALfile)
col.queen <- read.gal(GALfile)
unlink(GALfile)
summary(diffnb(col.queen, col.gal.nb))
}
\keyword{spatial}
