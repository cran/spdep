% Copyright 2002-3 by Roger S. Bivand
\name{LR.sarlm}
\alias{LR.sarlm}
\alias{LR1.sarlm}
\alias{Wald1.sarlm}
\alias{logLik.sarlm}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Likelihood ratio test}
\description{
  The \code{LR.sarlm()} function provides a likelihood ratio test for objects for which a \code{logLik()} function exists for their class, or for objects of class \code{logLik}. \code{LR1.sarlm()} and \code{Wald1.sarlm()} are used internally in \code{summary.sarlm()}, but may be accessed directly; they report the values respectively of LR and Wald tests for the absence of spatial dependence in spatial lag or error models.
}
\usage{
LR.sarlm(x, y)
logLik.sarlm(object, ...)
LR1.sarlm(object)
Wald1.sarlm(object)
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
  \code{LR1.sarlm} and \code{Wald1.sarlm} returm objects of class \code{htest}
}

\note{The numbers of degrees of freedom returned by \code{logLik.sarlm()} include nuisance parameters, that is the number of regression coefficients, plus sigma, plus spatial parameter esitmate(s).}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{logLik.lm}}, \code{\link{anova.sarlm}}}

\examples{
data(columbus)
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb),
  type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
LR.sarlm(mixed, error)
}
\keyword{spatial}
