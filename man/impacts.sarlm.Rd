\name{impacts.sarlm}
\Rdversion{1.1}
\alias{impacts.sarlm}
\alias{plot.sarlmImpact}
\alias{print.sarlmImpact}
\alias{summary.sarlmImpact}
\alias{print.summary.sarlmImpact}
\alias{HPDinterval.sarlmImpact}
\title{Impacts in spatial lag models}
\description{The calculation of impacts for spatial lag and spatial Durbin models is needed in order to interpret the regression coefficients correctly, because of the spillovers between the terms in these data generation processes (unlike the spatial error model).}
\usage{
impacts.sarlm(obj, \dots, tr, R = NULL, listw = NULL, useHESS = NULL, tol = 1e-06, empirical = FALSE)
\method{plot}{sarlmImpact}(x, ..., choice="direct", trace=FALSE, density=TRUE)
\method{print}{sarlmImpact}(x, ...)
\method{summary}{sarlmImpact}(object, ..., zstats=FALSE, short=FALSE)
\method{HPDinterval}{sarlmImpact}(obj, prob = 0.95, ..., choice="direct")
}
\arguments{
  \item{obj}{A sarlm spatial regression object created by \code{lagsarlm}; in \code{HPDinterval.sarlmImpact}, a sarlmImpact object}
  \item{\dots}{Arguments passed through to methods in the \pkg{coda} package}
  \item{tr}{A vector of traces of powers of the spatial weights matrix created using \code{trW}, for approximate impact measures; if not given, \code{listw} must be given for exact measures (for small to moderate spatial weights matrices); the traces must be for the same spatial weights as were used in fitting the spatial regression}
  \item{R}{If given, simulations are used to compute distributions for the impact measures, returned as \code{mcmc} objects}
  \item{listw}{If \code{tr} is not given, a spatial weights object as created by \code{nb2listw}; they must be the same spatial weights as were used in fitting the spatial regression}
  \item{useHESS}{Use the Hessian approximation (if available) even if the asymptotic coefficient covariance matrix is available; used for comparing methods}
  \item{tol}{Argument passed to \code{mvrnorm}: tolerance (relative to largest variance) for numerical lack of positive-definiteness in the coefficient covariance matrix}
  \item{empirical}{Argument passed to \code{mvrnorm} (default FALSE): if true, the coefficients and their covariance matrix specify the empirical not population mean and covariance matrix}
  \item{x, object}{sarlmImpact objects created by \code{impacts.sarlm}}
  \item{zstats}{default FALSE, if TRUE, also return z-values and p-values for the impacts based on the simulations}
  \item{short}{default FALSE, if TRUE passed to the print summary method to omit printing of the mcmc summaries}
  \item{choice}{One of three impacts: direct, indirect, or total}
  \item{trace}{Argument passed to \code{plot.mcmc}: plot trace plots}
  \item{density}{Argument passed to \code{plot.mcmc}: plot density plots}
  \item{prob}{Argument passed to \code{HPDinterval.mcmc}: a numeric scalar in the interval (0,1) giving the target probability content of the intervals}
}

\details{If called without \code{R} being set, the method returns the direct, indirect and total impacts for the variables in the model, for the variables themselves in tha spatial lag model case, for the variables and their spatial lags in the spatial Durbin (mixed) model case. The spatial lag impact measures are computed using eq. 2.46 (LeSage and Pace, 2009, p. 38), either using the exact dense matrix (when \code{listw} is given), or traces of powers of the weights matrix (when \code{tr} is given). When the traces are created by powering sparse matrices, the exact and the trace methods should give very similar results, unless the number of powers used is very small.

If \code{R} is given, simulations will be used to create distributions for the impact measures, provided that the fitted model object contains a coefficient covariance matrix. At present, this is only the case for the default eigenvalue method in \code{lagsarlm}, which provides an analytical covariance matrix, but alternatives will be provided for sparse matrix methods in the near future.

The simulations are stored as \code{mcmc} objects as defined in the \pkg{coda} package. The simulated values of the coefficients are checked to see that the spatial coefficient remains within its valid interval --- draws outside the interval are discarded.}

\value{
An object of class sarlmImpact.

If no simulation is carried out, the object returned is a list with:
\item{direct}{numeric vector}
\item{indirect}{numeric vector}
\item{total}{numeric vector}

If simulation is carried out, the object returned is a list with:

\item{res}{a list with three components as for the non-simulation case}
\item{sres}{a list with three \code{mcmc} matrices, for the direct, indirect and total impacts}
}

\references{LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton, pp. 33--42, 114--115; LeSage J and MM Fischer (2008) Spatial growth regressions: model specification, estimation and interpretation. Spatial Economic Analysis 3 (3), pp. 275--304.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{
\code{\link{trW}}, \code{\link{lagsarlm}}, \code{\link{nb2listw}}, \code{\link[MASS]{mvrnorm}}, \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{summary.mcmc}}, \code{\link[coda]{HPDinterval}}
}
\examples{
example(columbus)
listw <- nb2listw(col.gal.nb)
lobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw)
summary(lobj)
mobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed")
summary(mobj)
W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")
impacts.sarlm(lobj, listw=listw)
impacts.sarlm(lobj, tr=trMatc)
impacts.sarlm(lobj, tr=trMC)
summary(impacts.sarlm(lobj, tr=trMatc, R=1000), zstats=TRUE)
impacts.sarlm(mobj, listw=listw)
impacts.sarlm(mobj, tr=trMatc)
impacts.sarlm(mobj, tr=trMC)
summary(impacts.sarlm(mobj, tr=trMatc, R=1000), zstats=TRUE)
mobj1 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
method="Matrix", fdHess=TRUE)
summary(mobj1)
summary(impacts.sarlm(mobj1, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
summary(impacts.sarlm(mobj, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
mobj2 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
method="Matrix", fdHess=TRUE, optimHess=TRUE)
summary(impacts.sarlm(mobj2, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
mobj3 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
method="spam", fdHess=TRUE)
summary(impacts.sarlm(mobj3, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
data(boston)
Wb <- as(as_dgRMatrix_listw(nb2listw(boston.soi)), "CsparseMatrix")
trMatb <- trW(Wb, type="mult")
gp2mMi <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, nb2listw(boston.soi), type="mixed", method="Matrix", 
fdHess=TRUE, trs=trMatb)
summary(gp2mMi)
summary(impacts.sarlm(gp2mMi, tr=trMatb, R=1000), zstats=TRUE, short=TRUE)
data(house)
form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
   rooms + log(TLA) + beds + syear)
lobj <- lagsarlm(form, house, nb2listw(LO_nb), method="Matrix",
 fdHess=TRUE, trs=trMat)
summary(lobj)
loobj <- impacts.sarlm(lobj, tr=trMat, R=1000)
summary(loobj, zstats=TRUE, short=TRUE)
mobj <- lagsarlm(form, house, nb2listw(LO_nb), type="mixed",
 method="Matrix", fdHess=TRUE, trs=trMat)
summary(mobj)
moobj <- impacts.sarlm(mobj, tr=trMat, R=1000)
summary(moobj, zstats=TRUE, short=TRUE)
}
\keyword{spatial}

