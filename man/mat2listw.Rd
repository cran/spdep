% Copyright 2003 by Roger S. Bivand
\name{mat2listw}
\alias{mat2listw}

\title{Convert a square spatial weights matrix to a weights list object}
\description{
  The function converts a square spatial weights matrix to a weights list 
object, optionally adding region IDs from the row names of the matrix, as a 
sequence of numbers 1:nrow(x), or as given as an argument.
}
\usage{
mat2listw(x, row.names = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A square non-negative matrix with no NAs representing spatial 
weights}
  \item{row.names}{row names to use for region IDs}
}

\value{
A \code{listw} object with the following members:
  \item{style}{"M", meaning matrix style, underlying style unknown}
  \item{neighbours}{the derived neighbours list}
  \item{weights}{the weights for the neighbours derived from the matrix}
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

 
\seealso{\code{\link{nb2listw}}, \code{\link{nb2mat}}}
\examples{
data(columbus)
col005 <- dnearneigh(coords, 0, 0.5, attr(col.gal.nb, "region.id"))
summary(col005)
col005.w.mat <- nb2mat(col005, zero.policy=TRUE)
col005.w.b <- mat2listw(col005.w.mat)
summary(col005.w.b$neighbours)
diffnb(col005, col005.w.b$neighbours)
}
\keyword{spatial}
