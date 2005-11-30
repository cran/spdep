% Copyright 2002-4 by Roger S. Bivand
\name{errorsarlm}
\alias{errorsarlm}
%\alias{sar.error.f}
%\alias{sar.error.f.s}
\title{Spatial simultaneous autoregressive error model estimation}
\description{
  Maximum likelihood estimation of spatial simultaneous autoregressive
error models of the form:

\deqn{y = X \beta + u, u = \lambda W u + \varepsilon}{y = X beta + u, u = lambda W u + e}

where $\lambda$ is found by \code{optimize()}, or by \code{optim()} using method "L-BFGS-B" if argument optim=TRUE, first, and $\beta$ and other parameters by generalized least squares subsequently.
}
\usage{
errorsarlm(formula, data=list(), listw, na.action=na.fail, method="eigen",
  quiet=TRUE, zero.policy=FALSE, interval = c(-1, 0.999), tol.solve=1.0e-10, 
  tol.opt=.Machine$double.eps^0.5, control, optim=FALSE)
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
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw}, and "SparseM" for strictly symmetric weights lists of styles "B", "C" and "U", or made symmetric by similarity (Ord, 1975, Appendix C) if possible for styles "W" and "S", using code from the SparseM package to calculate the determinant.
}
  \item{quiet}{default=TRUE; if FALSE, reports function values during optimization.}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{errorsarlm()} to terminate with an error}
  \item{interval}{search interval for autoregressive parameter when not using method="eigen"; default is c(-1,1)}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=1.0e-10). This may be used if necessary to extract coefficient standard errors (for instance lowering to 1e-12), but errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{tol.opt}{the desired accuracy of the optimization - passed to \code{optim()} (default=square root of double precision machine tolerance, a larger root may be used if the warning: ERROR: ABNORMAL\_TERMINATION\_IN\_LNSRCH is seen, see help(boston) for an example)}
  \item{control}{A list of control parameters passed to \code{optim}, se details in \code{\link[stats]{optim}}}
  \item{optim}{If TRUE use experimental \code{optim} branch and control argument}

}
\details{
  The asymptotic standard error of $\lambda$ is only computed when
method=eigen, because the full matrix operations involved would be costly
for large n typically associated with the choice of method="SparseM". The
same applies to the coefficient covariance matrix. Taken as the
asymptotic matrix from the literature, it is typically badly scaled, being
block-diagonal, and with the elements involving lambda being very small,
while other parts of the matrix can be very large (often many orders
of magnitude in difference). It often happens that the \code{tol.solve}
argument needs to be set to a smaller value than the default, or the RHS variables can be centred or reduced in range.


}
\value{
  A list object of class \code{sarlm}
  \item{type}{"error"}
  \item{lambda}{simultaneous autoregressive error coefficient}
  \item{coefficients}{GLS coefficient estimates}
  \item{rest.se}{GLS coefficient standard errors (are equal to asymptotic
standard errors)}
  \item{LL}{log likelihood value at computed optimum}
  \item{s2}{GLS residual variance}
  \item{SSE}{sum of squared GLS errors}
  \item{parameters}{number of parameters estimated}
  \item{lm.model}{the \code{lm} object returned when estimating for $\lambda=0$}
  \item{method}{the method used to calculate the Jacobian}
  \item{call}{the call used to create this object}
  \item{residuals}{GLS residuals}
  \item{lm.target}{the \code{lm} object returned for the GLS fit}
  \item{fitted.values}{Difference between residuals and response variable}
  \item{ase}{TRUE if method=eigen}
  \item{formula}{model formula}
  \item{se.fit}{Not used yet}
  \item{lambda.se}{if ase=TRUE, the asymptotic standard error of $\lambda$}
  \item{LMtest}{NULL for this model}
  \item{zero.policy}{zero.policy for this model}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
  \item{aliased}{if not NULL, details of aliased variables}

The internal sar.error.* functions return the value of the log likelihood function at $\lambda$.
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
pp. 237-289.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{lm}}, \code{\link{lagsarlm}}, 
\code{\link{eigenw}}, \code{\link{asMatrixCsrListw}}, \code{\link{similar.listw}}, \code{\link{predict.sarlm}},
\code{\link{residuals.sarlm}}
}

\examples{
data(oldcol)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", quiet=FALSE)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", quiet=FALSE)
COL.errB.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"), method="eigen", quiet=FALSE)
COL.errW.SM <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="SparseM", quiet=FALSE)
summary(COL.errW.eig, correlation=TRUE)
summary(COL.errB.eig, correlation=TRUE)
summary(COL.errW.SM, correlation=TRUE)
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.err.NA <- errorsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 nb2listw(COL.nb), na.action=na.exclude)
COL.err.NA$na.action
COL.err.NA
resid(COL.err.NA)
}
\keyword{spatial}
