% Copyright 2002 by Roger S. Bivand
\name{LR.sarlm}
\alias{LR.sarlm}
\alias{logLik.sarlm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Likelihood ratio test}
\description{
  The function provides a likelihood ratio test for objects for which a \code{logLik()} function exists for their class, or for objects of class \code{logLik}.
}
\usage{
LR.sarlm(x, y)
logLik.sarlm(object, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{a \code{logLik} object or an object for which a \code{logLik()} function exists}
  \item{y}{a \code{logLik} object or an object for which a \code{logLik()} function exists}
  \item{object}{a \code{sarlm} object from \code{lagsarlm()} or \code{errorsarlm()}}
  \item{...}{other arguments to logLik()}
}
\value{
  \code{LR.sarlm()} returns an object of class \code{htest} with:
  \item{statistic}{value of statistic}
  \item{parameter}{degrees of freedom}
  \item{p.value}{Probability value}
  \item{estimate}{Log likelihood values of compared models}
  \item{method}{"Likelihood ratio for spatial linear models"}

  \code{logLik.sarlm()} returns an object of class \code{logLik}
}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{logLik.lm}}}

\examples{
data(oldcol)
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=COL.OLD, nb2listw(COL.nb),
  type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=COL.OLD, nb2listw(COL.nb))
LR.sarlm(mixed, error)
}
\keyword{spatial}
