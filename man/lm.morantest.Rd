% Copyright 2001-5 by Roger S. Bivand
\name{lm.morantest}
\alias{lm.morantest}
\alias{listw2U}
\title{Moran's I test for residual spatial autocorrelation}
\description{
  Moran's I test for spatial autocorrelation in residuals from an estimated linear model (\code{lm()}). The helper function \code{listw2U()} constructs a weights list object corresponding to the sparse matrix \eqn{\frac{1}{2} ( \mathbf{W} + \mathbf{W}'}{1/2 (W + W')}
}
\usage{
lm.morantest(model, listw, zero.policy=FALSE, 
	    alternative = "greater", spChk=NULL, resfun=weighted.residuals)
listw2U(listw)
}
\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm}; weights
may be specified in the \code{lm} fit, but offsets should not be used}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without
neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis,
must be one of "greater" (default), "less" or "two.sided".}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{resfun}{default: weighted.residuals; the function to be used to extract residuals from the \code{lm} object, may be \code{residuals}, \code{weighted.residuals}, \code{rstandard}, or \code{rstudent}}
}
\value{
A list with class \code{htest} containing the following components:
  \item{statistic}{the value of the standard deviate of Moran's I.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed Moran's I, its expectation and
variance under the method assumption.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 203, }
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{lm.LMtests}}, \code{\link{lm}}}

\note{The examples also show how to use permutation bootstrap if desired.}

\examples{
data(oldcol)
oldcrime1.lm <- lm(CRIME ~ 1, data = COL.OLD)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="W"))
lm.LMtests(oldcrime.lm, nb2listw(COL.nb, style="W"))
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="S"))
lm.morantest(oldcrime1.lm, nb2listw(COL.nb, style="W"))
moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
oldcrime.wlm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD, weights = I(1/AREA))
lm.morantest(oldcrime.wlm, nb2listw(COL.nb, style="W"), resfun=weighted.residuals)
lm.morantest(oldcrime.wlm, nb2listw(COL.nb, style="W"), resfun=rstudent)
if (require(boot)) {
  oldcrime.lmx <- lm(CRIME ~ HOVAL + INC, data = COL.OLD, x=TRUE)
  listw <- nb2listw(COL.nb, style="W")
  MoraneI.boot <- function(var, i, ...) {
    var <- var[i]
    lmres <- lm(var ~ oldcrime.lmx$x - 1)
    return(moran(x=residuals(lmres), ...)$I)
  }
  boot1 <- boot(residuals(oldcrime.lmx), statistic=MoraneI.boot, R=499,
    sim="permutation", listw=listw, n=length(listw$neighbours), S0=Szero(listw))
  zi <- (boot1$t0 - mean(boot1$t))/sqrt(var(boot1$t))
  boot1
  plot(boot1)
  cat("Bootstrap permutation standard deviate:", zi, "\n\n")
  lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="W"))
}
}
\keyword{spatial}
