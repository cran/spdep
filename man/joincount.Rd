% Copyright 2001 by Roger S. Bivand
\name{joincount}
\alias{joincount}
\title{Count joins between same colours}
\description{
  A simple function to compute the generalised BB join count statistic for k-coloured factors, called by \code{joincount.test} and \code{joincount.mc};
\deqn{BB = \sum_{i=1}^{n}\sum_{j=1}^{n}w_{ij} x_i x_j}
}
\usage{
joincount(dums, listw)
}

\arguments{
  \item{dums}{a matrix of binary dummy variables created from a factor}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
}

\value{
a vector of the raw weighted count sums for each level of the factor
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 17.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{joincount.test}}, \code{\link{joincount.mc}}}

\keyword{spatial}
