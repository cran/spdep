% Copyright 2001 by Roger S. Bivand
\name{lag.listw}
\alias{lag.listw}
\title{Spatial lag of a numeric vector}
\description{
 Using a \code{listw} sparse representation of a spatial weights matrix, compute the lag vector \eqn{V x}
}
\usage{
lag.listw(listw, x, zero.policy=FALSE)
}

\arguments{
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
}
\value{
  a numeric vector the same length as x
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link[spweights]{nb2listw}}}

\examples{
data(oldcol)
Vx <- lag.listw(nb2listw(COL.nb, style="W"), COL.OLD$CRIME)
plot(Vx, COL.OLD$CRIME)
require(stepfun)
plot(ecdf(COL.OLD$CRIME))
plot(ecdf(Vx), add=TRUE, col.points="red", col.hor="red")
}
\keyword{spatial}
