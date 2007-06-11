% Copyright 2001-4 by Roger S. Bivand
\encoding{latin1}
\name{moran.mc}
\alias{moran.mc}
\title{Permutation test for Moran's I statistic}
\description{
 A permutation test for Moran's I statistic calculated by using nsim random permutations of x for the given spatial weighting scheme, to establish the rank of the observed statistic in relation to the nsim simulated values. The examples show how \code{boot(sim="permutation")} can replicate this function (thanks to Virgilio Gómez Rubio and the DCluster package).
}
\usage{
moran.mc(x, listw, nsim, zero.policy=FALSE, alternative="greater",
 na.action=na.fail, spChk=NULL)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of "greater" (default), or "less".}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted. \code{na.pass} is not permitted because it is meaningless in a permutation test.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}}
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
colw <- nb2listw(COL.nb, style="W")
nsim <- 99
set.seed(1234)
sim1 <- moran.mc(COL.OLD$CRIME, listw=colw, nsim=nsim)
sim1
mean(sim1$res[1:nsim])
var(sim1$res[1:nsim])
summary(sim1$res[1:nsim])
MoranI.boot <- function(var, i, ...) {
var <- var[i]
return(moran(x=var, ...)$I)
}
set.seed(1234)
library(boot)
boot1 <- boot(COL.OLD$CRIME, statistic=MoranI.boot, R=nsim,
 sim="permutation", listw=colw, n=nrow(COL.OLD), S0=Szero(colw))
boot1
plot(boot1)
mean(boot1$t)
var(boot1$t)
summary(boot1$t)
colold.lags <- nblag(COL.nb, 3)
set.seed(1234)
sim2 <- moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"), nsim=nsim)
summary(sim2$res[1:nsim])
sim3 <- moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), nsim=nsim)
summary(sim3$res[1:nsim])
}
\keyword{spatial}
