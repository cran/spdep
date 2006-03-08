\name{GMerrorsar}
\alias{GMerrorsar}
\alias{residuals.gmsar}
\alias{deviance.gmsar}
\alias{coef.gmsar}
\alias{fitted.gmsar}
\alias{print.gmsar}
\alias{summary.gmsar}
\alias{LR1.gmsar}
\alias{logLik.gmsar}
\alias{print.summary.gmsar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial simultaneous autoregressive error model estimation by GMM}
\description{
  An implementation of Kelejian and Prucha's generalised moments estimator for the autoregressive parameter in a spatial model.
}
\usage{
GMerrorsar(formula, data = list(), listw, na.action = na.fail, zero.policy = FALSE, return_LL = TRUE, control = list(), verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{GMerrorsar()} to terminate with an error}
  \item{control}{A list of control parameters. See details in \code{\link[stats]{optim}}.}
  \item{return_LL}{default TRUE, if FALSE, do not try to calculate the log likelihood of the function for the fitted model values --- see details}
  \item{verbose}{default=FALSE; if TRUE, reports function values during optimization.}
}
\details{
  When the control list is set with care, the function will converge to values close to the ML estimator without requiring computation of the Jacobian, the most resource-intensive part of ML estimation. For moderately sized data sets with hundreds of observations, but not many thousands, the Jacobian is computed once to give the likelihood of the fitted model, allowing a test against the model with no spatial dependence.

Note that the fitted() function for the output object assumes that the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564). This differs from other software, including GeoDa, which does not use knowledge of the response 
variable in making predictions for the fitting data.

}
\value{
  A list object of class \code{gmsar}
  \item{lambda}{simultaneous autoregressive error coefficient}
  \item{coefficients}{GMM coefficient estimates}
  \item{rest.se}{GMM coefficient standard errors}
  \item{s2}{GMM residual variance}
  \item{SSE}{sum of squared GMM errors}
  \item{parameters}{number of parameters estimated}
  \item{lm.model}{the \code{lm} object returned when estimating for $\lambda=0$}
  \item{call}{the call used to create this object}
  \item{residuals}{GMM residuals}
  \item{lm.target}{the \code{lm} object returned for the GMM fit}
  \item{fitted.values}{Difference between residuals and response variable}
  \item{formula}{model formula}
  \item{aliased}{if not NULL, details of aliased variables}
  \item{zero.policy}{zero.policy for this model}
  \item{LL}{log likelihood value at computed optimum}
  \item{na.action}{(possibly) named vector of excluded or omitted observations if non-default na.action argument used}
}
\references{Kelejian, H. H., and Prucha, I. R., 1999. A Generalized Moments Estimator for the Autoregressive Parameter in a Spatial Model. International Economic Review, 40, pp. 509--533; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York.}
\author{Luc Anselin and Roger Bivand}

\seealso{\code{\link[stats]{optim}}, \code{\link{errorsarlm}}}
\examples{
data(oldcol)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb, style="W"), method="eigen")
summary(COL.errW.eig)
COL.errW.GM <- GMerrorsar(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb, style="W"))
summary(COL.errW.GM)
data(NY_data)
esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata, listw=listw_NY, family="SAR", method="full")
summary(esar1f)
esar1gm <- GMerrorsar(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata, listw=listw_NY)
summary(esar1gm)
}
\keyword{spatial}

