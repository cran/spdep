% Copyright 2001-2 by Roger S. Bivand
\name{moran.mc}
\alias{moran.mc}
\title{Permutation test for Moran's I statistic}
\description{
 A permutation test for Moran's I statistic calculated by using nsim random permutations of x for the given spatial weighting scheme, to establish the rank of the observed statistic in relation to the nsim simulated values.
}
\usage{
moran.mc(x, listw, nsim, zero.policy=FALSE, alternative="greater")
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.}
}

\value{
A list with class \code{htest} and \code{mc.sim} containing the following components:
  \item{statistic}{the value of the observed Moran's I.}
  \item{parameter}{the rank of the observed Moran's I.}
  \item{p.value}{the pseudo p-value of the test.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data, and the number of simulations.}
  \item{res}{nsim simulated values of statistic, final value is observed statistic}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 63-5.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{moran}}, \code{\link{moran.test}}}

\examples{
data(oldcol)
sim1 <- moran.mc(COL.OLD$CRIME, nb2listw(COL.nb, style="W"), nsim=99)
mean(sim1$res)
var(sim1$res)
summary(sim1$res)
colold.lags <- nblag(COL.nb, 3)
sim2 <- moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]], style="W"),
 nsim=99)
summary(sim2$res)
sim3 <- moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]], style="W"),
 nsim=99)
summary(sim3$res)
}
\keyword{spatial}
