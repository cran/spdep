% Copyright 2001 by Roger S. Bivand
\name{is.symmetric.nb}
\alias{is.symmetric.nb}
\alias{sym.attr.nb}
\alias{make.sym.nb}

\title{Test a neighbours list for symmetry}
\description{
Checks a neighbours list for symmetry/transitivity (if i is a neighbour of j,
then j is a neighbour of i). This holds for distance and contiguity based
neighbours, but not for k-nearest neighbours. The helper function
\code{sym.attr.nb()} calls \code{is.symmetric.nb()} to set the \code{sym}
attribute if needed, and \code{make.sym.nb} makes a non-symmetric list symmetric by adding neighbors.
}
\usage{
is.symmetric.nb(nb, verbose = TRUE, force = FALSE)
sym.attr.nb(nb)
make.sym.nb(nb)
}

\arguments{
  \item{nb}{an object of class \code{nb} with a list of integer vectors containing neighbour region number ids.}
  \item{verbose}{if TRUE prints non-matching pairs}
  \item{force}{do not respect a neighbours list \code{sym} attribute and test anyway}
}
\value{
TRUE if symmetric, FALSE if not
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{read.gal}}}

\examples{
data(columbus)
print(is.symmetric.nb(col.gal.nb, verbose=TRUE, force=TRUE))
k4 <- knn2nb(knearneigh(coords, k=4), row.names=rownames(columbus))
k4 <- sym.attr.nb(k4)
print(is.symmetric.nb(k4))
k4.sym <- make.sym.nb(k4)
print(is.symmetric.nb(k4.sym))
}
\keyword{spatial}
