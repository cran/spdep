% Copyright 2001 by Roger S. Bivand
\name{card}
\alias{card}
\title{Cardinalities for neighbours lists}
\description{
The function tallies the numbers of neighbours of regions in the neighbours
list.
}
\usage{
card(nb)
}
\arguments{
  \item{nb}{a neighbours list object of class \code{nb}}
}
\value{
An integer vector of the numbers of neighbours of regions in the neighbours
list.
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}}

\examples{
example(columbus)
table(card(col.gal.nb))
}
\keyword{spatial}

