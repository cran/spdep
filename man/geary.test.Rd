% Copyright 2001 by Roger S. Bivand
\name{geary.test}
\alias{geary.test}
\title{Geary's C test for spatial autocorrelation}
\description{
   Geary's test for spatial autocorrelation using a spatial weights matrix in weights list form. The assumptions underlying the test are sensitive to the form of the graph of neighbour relationships and other factors, and results may be checked against those of \code{geary.mc} permutations.
}
\usage{
geary.test(x, listw, randomisation=TRUE, zero.policy=FALSE,
    alternative="less", spChk=NULL)
}

\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{randomisation}{variance of I calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of less (default), greater or two.sided.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
}

\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Geary's C.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed Geary's C, its expectation and 
 variance under the method assumption.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the assumption used for calculating 
 the standard deviate.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{Probability values one-tailed, as the hypothesis most often used is no 
autocorrelation against positive autocorrelation.}

\seealso{\code{\link{geary}}, \code{\link{geary.mc}}}

\examples{
data(oldcol)
geary.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb, style="W"))
geary.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
colold.lags <- nblag(COL.nb, 3)
geary.test(spNamedVec("CRIME", COL.OLD), nb2listw(colold.lags[[2]],
 style="W"))
geary.test(spNamedVec("CRIME", COL.OLD), nb2listw(colold.lags[[3]],
 style="W"), alternative="greater")
}
\keyword{spatial}
