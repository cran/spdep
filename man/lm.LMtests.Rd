% Copyright 2001-2 by Roger S. Bivand
\name{lm.LMtests}
\alias{lm.LMtests}
\alias{print.LMtestlist}
\alias{tracew}
\title{Lagrange Multiplier diagnostics for spatial dependence in linear models}
\description{
The function reports the estimates of tests chosen among five statistics for
testing for spatial dependence in linear models. The statistics are
the simple LM test for error dependence (LMerr), the simple LM test
for a missing spatially lagged dependent variable (LMlag), variants
of these robust to the presence of the other (RLMerr, RLMlag - RLMerr
tests for error dependence in the possible presence of a missing lagged
dependent variable, RLMlag the other way round), and a portmanteau test
(SARMA, in fact LMerr + RLMlag).
}
\usage{
lm.LMtests(model, listw, zero.policy=FALSE, test="LMerr", spChk=NULL)
print.LMtestlist(x, ...)
tracew(listw)
}

\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm}; weights
and offsets should not be used}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw},
expected to be row-standardised (W-style)}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without
neighbours, if FALSE assign NA}
  \item{test}{a character vector of tests requested chosen from LMerr, LMlag,
RLMerr, RLMlag, SARMA; test="all" computes all the tests.}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{x}{object to be printed}
  \item{...}{printing arguments to be passed through}
}

\details{
The two types of dependence are for spatial lag $\rho$ and spatial error
$\lambda$:

\deqn{
\mathbf{y} = \mathbf{X \beta} + \rho \mathbf{W_{(1)} y} + \mathbf{u},
}{y = X beta + rho W1 y + u}
\deqn{
\mathbf{u} = \lambda \mathbf{W_{(2)} u} + \mathbf{e}
}{u = lambda W2 u + e}

where \eqn{\mathbf{e}}{e} is a well-behaved, uncorrelated error
term. Tests for a missing spatially lagged dependent variable test
that \eqn{\rho = 0}{rho = 0}, tests for spatial autocorrelation of
the error \eqn{\mathbf{u}}{u} test whether \eqn{\lambda = 0}{lambda =
0}. \eqn{\mathbf{W}}{W} is a spatial weights matrix; for the tests used
here they are identical.
}

\value{
A list of class \code{LMtestlist} of \code{htest} objects, each with:
  \item{statistic}{the value of the Lagrange Multiplier test.}
  \item{parameter}{number of degrees of freedom}
  \item{p.value}{the p-value of the test.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Anselin, L. 1988 Spatial econometrics: methods and
models. (Dordrecht: Kluwer); Anselin, L., Bera, A. K., Florax, R. and
Yoon, M. J. 1996 Simple diagnostic tests for spatial dependence. Regional
Science and Urban Economics, 26, 77--104.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Andrew Bernat}

\seealso{\code{\link{lm}}}

\examples{
data(oldcol)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
summary(oldcrime.lm)
lm.LMtests(oldcrime.lm, nb2listw(COL.nb), test=c("LMerr", "LMlag", "RLMerr",
  "RLMlag", "SARMA"))
}
\keyword{spatial}
