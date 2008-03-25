\name{spautolm}
\alias{spautolm}
\alias{residuals.spautolm}
\alias{deviance.spautolm}
\alias{coef.spautolm}
\alias{fitted.spautolm}
\alias{print.spautolm}
\alias{summary.spautolm}
\alias{LR1.spautolm}
\alias{logLik.spautolm}
\alias{print.summary.spautolm}
%\alias{CAR}
%\alias{SAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial conditional and simultaneous autoregression model estimation}
\description{
 Function taking family and weights arguments for spatial autoregression model estimation by Maximum Likelihood, using full matrix methods, not suited to large data sets with thousands of observations. With one of the sparse matrix methods, larger numbers of observations can be handled, but the \code{interval=} argument should be set. The implementation is GLS using the single spatial coefficient value, here termed lambda, found by line search using \code{optimize} to maximise the log likelihood.
}
\usage{
spautolm(formula, data = list(), listw, weights,
 na.action=na.fail, verbose = FALSE, tol.opt = .Machine$double.eps^(2/3),
 family = "SAR", method="full", interval=c(-1,0.999),
 zero.policy = FALSE, super=FALSE, Matrix_intern=TRUE, %cholAlloc=NULL, 
 tol.solve=.Machine$double.eps, llprof=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{weights}{an optional vector of weights to be used in the fitting process}
  \item{verbose}{default=FALSE; if TRUE, reports function values during optimization.}
  \item{tol.opt}{the desired accuracy of the optimization - passed to \code{optimize()}}
  \item{family}{character string: either \code{"SAR"} or \code{"CAR"} for simultaneous or conditional autoregressions; \code{"SMA"} for spatial moving average added thanks to Jielai Ma - \code{"SMA"} is only implemented for method=\code{"full"} because it necessarily involves dense matrices}
  \item{method}{character string: default \code{"full"} for use of full matrices, \code{"Matrix"} or \code{"spam"} for sparse matrices (restricted to spatial weights symmetric or similar to symmetric) using methods in either the Matrix or spam packages. Note that users should change from spam to Matrix}
  \item{interval}{search interval for autoregressive parameter when not using method="eigen"; default is c(-1,0.999), \code{optimize} will reset NA/NaN to a bound and gives a warning when the interval is poorly set; method="Matrix" will attempt to search for an appropriate interval}
  \item{zero.policy}{Include list of no-neighbour observations in output if TRUE --- otherwise zero.policy is handled within the listw argument}
  \item{super}{logical argument to \code{\link[Matrix]{Cholesky}}, only used for method="Matrix"; at present only super=FALSE is supported. It should only be used for positive lambda. If NULL, the Cholesky decompisition of the matrix is undertaken at each step.}
  \item{Matrix_intern}{default TRUE, if FALSE, does not use vectorised Matrix internal function.}
%  \item{cholAlloc}{control arguments for memory allocation in \code{chol} for
% sparse matrices in method="SparseM": if NULL, default values are used
% (larger than SparseM defaults), otherwise a list with elements: nsubmax,
% nnzlmax, and tmpmax, for example cholAlloc=list(nsubmax=25000,
% nnzlmax=100000, tmpmax=25000)}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=double precision machine tolerance). Errors in \code{solve()} may constitute indications of poorly scaled variables: if the variables have scales differing much from the autoregressive coefficient, the values in this matrix may be very different in scale, and inverting such a matrix is analytically possible by definition, but numerically unstable; rescaling the RHS variables alleviates this better than setting tol.solve to a very small value}
  \item{llprof}{default NULL, can either be an integer for method="full" only, to divide the feasible range into llprof points, or a sequence of values, at which to evaluate the likelihood function}
}
\details{
  This implementation is based on \code{\link[MASS]{lm.gls}} and \code{\link{errorsarlm}}. In particular, the function does not (yet) prevent asymmetric spatial weights being used with "CAR" family models. It appears that both numerical issues (convergence in particular) and uncertainties about the exact spatial weights matrix used make it difficult to reproduce Cressie and Chan's 1989 results, also given in Cressie 1993.

Note that the fitted() function for the output object assumes that the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564). This differs from other software, including GeoDa, which does not use knowledge of the response 
variable in making predictions for the fitting data.

}
\value{
  A list object of class \code{spautolm}:
  \item{fit}{a list, with items:
    \item{coefficients}{ML coefficient estimates}
    \item{SSE}{ML sum of squared errors}
    \item{s2}{ML residual variance}
    \item{imat}{ML coefficient covariance matrix}
    \item{signal\_trend}{non-spatial component of fitted.values}
    \item{signal\_stochastic}{spatial component of fitted.values}
    \item{fitted.values}{sum of non-spatial and spatial components of fitted.values}
    \item{residuals}{difference between observed and fitted values}
  }
  \item{lambda}{ML autoregressive coefficient}
  \item{LL}{log likelihood for fitted model}
  \item{LL0}{log likelihood for model with lambda=0}
  \item{call}{the call used to create this object}
  \item{parameters}{number of parameters estimated}
  \item{aliased}{if not NULL, details of aliased variables}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
  \item{llprof}{if not NULL, a list with components lambda and ll of equal length}
}
\references{Cliff, A. D., Ord, J. K. 1981 \emph{Spatial processes}, Pion;
Ord, J. K. 1975 Estimation methods for models of spatial interaction,
\emph{Journal of the American Statistical Association}, 70, 120-126; Waller, L. A., Gotway, C. A. 2004 \emph{Applied spatial statistics for public health}, Wiley, Hoboken, NJ, 325-380; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York, 548-568; Ripley, B. D. 1981 \emph{Spatial statistics}, Wiley, New York, 88-95.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{The standard errors given in Waller and Gotway (2004) are adjusted for the numbers of parameters estimated, and may be reproduced by using the additional argument \code{adj.se=TRUE} in the \code{summary} method. In addition, the function returns fitted values and residuals as given by Cressie (1993) p. 564.}

\seealso{\code{\link[stats]{optimize}}, \code{\link{errorsarlm}}}
\examples{
\dontrun{
esarI <- errorsarlm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
 zero.policy=TRUE)
summary(esarI)
esarIa <- spautolm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
 family="SAR")
summary(esarIa)
esarIV <- errorsarlm(ft.SID74 ~ ft.NWBIR74, data=mdata, listw=sids.nhbr.listw,
 zero.policy=TRUE)
summary(esarIV)
esarIVa <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata, listw=sids.nhbr.listw,
 family="SAR")
summary(esarIVa)
esarIaw <- spautolm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
 weights=BIR74, family="SAR")
summary(esarIaw)
esarIIaw <- spautolm(ft.SID74 ~ both - 1, data=mdata, listw=sids.nhbr.listw,
 weights=BIR74, family="SAR")
summary(esarIIaw)
esarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata,
 listw=sids.nhbr.listw, weights=BIR74, family="SAR")
summary(esarIVaw)
ecarIaw <- spautolm(ft.SID74 ~ 1, data=mdata.4, listw=sids.nhbr.listw.4,
 weights=BIR74, family="CAR")
summary(ecarIaw)
ecarIIaw <- spautolm(ft.SID74 ~ both - 1, data=mdata.4,
 listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIIaw)
ecarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata.4,
 listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
summary(ecarIVaw)
nc_LL$fitIV <- append(fitted.values(ecarIVaw), NA, outl-1)
spplot(nc_LL, c("fitIV"), cuts=12) # Cressie 1993, p. 565
}
example(NY_data)
lm0 <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata)
summary(lm0)
lm0w <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata, weights=POP8)
summary(lm0w)
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(esar0)
system.time(esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="full", verbose=TRUE))
summary(esar1f)
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="Matrix", verbose=TRUE))
summary(esar1M)
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="Matrix", verbose=TRUE,
 super=NULL))
summary(esar1M)
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="Matrix", verbose=TRUE,
 Matrix_intern=FALSE))
summary(esar1M)
system.time(esar1s <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="spam", verbose=TRUE))
summary(esar1s)
esar1wf <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="full")
summary(esar1wf)
system.time(esar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, weights=POP8, family="SAR", method="Matrix"))
summary(esar1wM)
system.time(esar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, weights=POP8, family="SAR", method="Matrix",
 super=FALSE))
summary(esar1wM)
esar1ws <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="spam")
summary(esar1ws)
ecar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="CAR", method="full")
summary(ecar1f)
system.time(ecar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="CAR", method="Matrix"))
summary(ecar1M)
system.time(ecar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="CAR", method="Matrix", super=FALSE))
summary(ecar1M)
ecar1s <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="CAR", method="spam")
summary(ecar1s)
ecar1wf <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=nydata$POP8, family="CAR", method="full")
summary(ecar1wf)
system.time(ecar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, weights=POP8, family="CAR", method="Matrix"))
summary(ecar1wM)
system.time(ecar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, weights=POP8, family="CAR", method="Matrix",
 super=FALSE))
summary(ecar1wM)
ecar1ws <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="CAR", method="spam")
summary(ecar1ws)
data(oldcol)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"))
summary(COL.errW.eig)
COL.errW.sar <- spautolm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"))
summary(COL.errW.sar)
data(boston)
gp1 <- spautolm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2)
 + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
 data=boston.c, nb2listw(boston.soi), family="SMA")
summary(gp1)

}
\keyword{spatial}

