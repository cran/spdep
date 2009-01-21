% Copyright 2002-8 by Roger S. Bivand
\name{lagsarlm}
\alias{lagsarlm}
\title{Spatial simultaneous autoregressive lag model estimation}
\description{
  Maximum likelihood estimation of spatial simultaneous autoregressive
lag and mixed models of the form:

\deqn{y = \rho W y + X \beta + \varepsilon}{y = rho W y + X beta + e}

where \eqn{\rho}{rho} is found by \code{optimize()} first, and \eqn{\beta}{beta} and other parameters by generalized least squares subsequently (one-dimensional search using optim performs badly on some platforms). In the mixed model, the spatially lagged independent variables are added to X.
}
\usage{
lagsarlm(formula, data=list(), listw, na.action=na.fail,
  type="lag", method="eigen", quiet=TRUE,
  zero.policy=FALSE, interval = c(-1, 0.999), tol.solve=1.0e-10, 
  tol.opt=.Machine$double.eps^0.5%, cholAlloc=NULL
)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{type}{default "lag", may be set to "mixed"; when "mixed", the lagged intercept is dropped for spatial weights style "W", that is row-standardised weights, but otherwise included}
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw}, and "spam" or "Matrix" for strictly symmetric weights lists of styles "B" and "C", or made symmetric by similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S", using code from the spam or Matrix packages to calculate the determinant.
}
  \item{quiet}{default=TRUE; if FALSE, reports function values during optimization.}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{lagsarlm()} to terminate with an error}
  \item{interval}{search interval for autoregressive parameter when not using method="eigen"; default is c(-1,1); method="Matrix" will attempt to search for an appropriate interval}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=1.0e-10). This may be used if necessary to extract coefficient standard errors (for instance lowering to 1e-12), but errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{tol.opt}{the desired accuracy of the optimization - passed to \code{optimize()} (default=square root of double precision machine tolerance)}
%  \item{cholAlloc}{control arguments for memory allocation in \code{chol}
% for sparse matrices in method="SparseM": if NULL, default values are
% used (larger than SparseM defaults), otherwise a list with elements:
% nsubmax, nnzlmax, and tmpmax, for example cholAlloc=list(nsubmax=25000,
% nnzlmax=100000, tmpmax=25000)}
}
\details{
  The asymptotic standard error of \eqn{\rho}{rho} is only computed when
method=eigen, because the full matrix operations involved would be costly
for large n typically associated with the choice of method="spam" or "Matrix". The
same applies to the coefficient covariance matrix. Taken as the
asymptotic matrix from the literature, it is typically badly scaled, and with the elements involving \eqn{\rho}{rho} being very small,
while other parts of the matrix can be very large (often many orders
of magnitude in difference). It often happens that the \code{tol.solve}
argument needs to be set to a smaller value than the default, or the RHS variables can be centred or reduced in range.

Note that the fitted() function for the output object assumes that the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564). This differs from other software, including GeoDa, which does not use knowledge of the response 
variable in making predictions for the fitting data.
}

\value{
  A list object of class \code{sarlm}
  \item{type}{"lag" or "mixed"}
  \item{rho}{simultaneous autoregressive lag coefficient}
  \item{coefficients}{GLS coefficient estimates}
  \item{rest.se}{asymptotic standard errors if ase=TRUE}
  \item{LL}{log likelihood value at computed optimum}
  \item{s2}{GLS residual variance}
  \item{SSE}{sum of squared GLS errors}
  \item{parameters}{number of parameters estimated}
  \item{lm.model}{the \code{lm} object returned when estimating for \eqn{\rho=0}{rho=0}}
  \item{method}{the method used to calculate the Jacobian}
  \item{call}{the call used to create this object}
  \item{residuals}{GLS residuals}
  \item{lm.target}{the \code{lm} object returned for the GLS fit}
  \item{fitted.values}{Difference between residuals and response variable}
  \item{se.fit}{Not used yet}
  \item{formula}{model formula}
  \item{ase}{TRUE if method=eigen}
  \item{LLs}{if ase=FALSE (for method="spam" or "Matrix"), the log likelihood values of
models estimated dropping each of the independent variables in turn, used
in the summary function as a substitute for variable coefficient
significance tests}
  \item{rho.se}{if ase=TRUE, the asymptotic standard error of \eqn{\rho}{rho}}
  \item{LMtest}{if ase=TRUE, the Lagrange Multiplier test for the absence
of spatial autocorrelation in the lag model residuals}
  \item{zero.policy}{zero.policy for this model}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}

The internal sar.lag.mixed.* functions return the value of the log likelihood function at \eqn{\rho}{rho}.

}
\references{Cliff, A. D., Ord, J. K. 1981 \emph{Spatial processes}, Pion;
Ord, J. K. 1975 Estimation methods for models of spatial interaction,
\emph{Journal of the American Statistical Association}, 70, 120-126;
Anselin, L. 1988 \emph{Spatial econometrics: methods and models.}
(Dordrecht: Kluwer); Anselin, L. 1995 SpaceStat, a software program for
the analysis of spatial data, version 1.80. Regional Research Institute,
West Virginia University, Morgantown, WV (\url{www.spacestat.com});
Anselin L, Bera AK (1998) Spatial dependence in linear regression models
with an introduction to spatial econometrics. In: Ullah A, Giles DEA
(eds) Handbook of applied economic statistics. Marcel Dekker, New York,
pp. 237-289; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}, with thanks to Andrew 
Bernat for contributions to the asymptotic standard error code.}

\seealso{\code{\link{lm}}, \code{\link{errorsarlm}}, 
\code{\link{eigenw}}, %\code{\link{asMatrixCsrListw}}, 
\code{\link{predict.sarlm}},
\code{\link{residuals.sarlm}}
}

\examples{
data(oldcol)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", quiet=FALSE)
summary(COL.lag.eig, correlation=TRUE)
system.time(COL.lag.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="Matrix", quiet=FALSE))
summary(COL.lag.M, correlation=TRUE)
system.time(COL.lag.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="spam", quiet=FALSE))
summary(COL.lag.M, correlation=TRUE)
COL.lag.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"))
summary(COL.lag.B, correlation=TRUE)
COL.mixed.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"), type="mixed", tol.solve=1e-9)
summary(COL.mixed.B, correlation=TRUE)
COL.mixed.W <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), type="mixed")
summary(COL.mixed.W, correlation=TRUE)
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.lag.NA <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 nb2listw(COL.nb), na.action=na.exclude, tol.opt=.Machine$double.eps^0.4)
COL.lag.NA$na.action
COL.lag.NA
resid(COL.lag.NA)
}
\keyword{spatial}
