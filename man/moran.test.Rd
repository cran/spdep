% Copyright 2001 by Roger S. Bivand
\name{moran.test}
\alias{moran.test}
\title{Moran's I test for spatial autocorrelation}
\description{
 Moran's test for spatial autocorrelation using a spatial weights matrix in weights list form. The assumptions underlying the test are sensitive to the form of the graph of neighbour relationships and other factors, and results may be checked against those of \code{moran.mc} permutations.
}
\usage{
moran.test(x, listw, randomisation=TRUE, zero.policy=FALSE,
	alternative="greater", rank = FALSE, spChk=NULL)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{randomisation}{variance of I calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.}
  \item{rank}{logical value - default FALSE for continuous variables, if TRUE, uses the adaptation of Moran's I for ranks suggested by Cliff and Ord (1981, p. 46)}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
}


\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Moran's I.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed Moran's I, its expectation and variance under the method assumption.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the assumption used for calculating the standard deviate.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\note{Var(I) is taken from Goodchild's CATMOG 47,
see also Upton & Fingleton (1985) p. 171; it agrees with SpaceStat,
see Tutorial workbook Chapter 22; VI is as given by Cliff and Ord minus the
square of EI.}

\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{moran}}, \code{\link{moran.mc}}}

\examples{
data(oldcol)
moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb, style="W"))
moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb, style="B"))
moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb, style="C"))
moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb, style="S"))
moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
colold.lags <- nblag(COL.nb, 3)
moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(colold.lags[[2]],
 style="W"))
moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(colold.lags[[3]],
 style="W"))
ranks <- rank(COL.OLD$CRIME)
names(ranks) <- rownames(COL.OLD)
moran.test(ranks, nb2listw(COL.nb, style="W"), rank=TRUE)
}
\keyword{spatial}
