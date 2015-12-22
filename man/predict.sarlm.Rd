% Copyright 2002 by Roger S. Bivand
\name{predict.sarlm}
\alias{predict.sarlm}
\alias{print.sarlm.pred}
\alias{as.data.frame.sarlm.pred}

\title{Prediction for spatial simultaneous autoregressive linear
model objects}
\description{
  \code{predict.sarlm()} calculates predictions as far as is at present
possible for for spatial simultaneous autoregressive linear
model objects, using Haining's terminology for decomposition into
trend, signal, and noise --- see reference.
}
\usage{
\method{predict}{sarlm}(object, newdata = NULL, listw = NULL, type = NULL, all.data = FALSE,
 zero.policy = NULL, legacy = TRUE, legacy.mixed = FALSE, power = NULL, order = 250,
 tol = .Machine$double.eps^(3/5), spChk = NULL, ...)
\method{print}{sarlm.pred}(x, ...)
\method{as.data.frame}{sarlm.pred}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{\code{sarlm} object returned by \code{lagsarlm}, 
\code{errorsarlm} or \code{sacsarlm}}
  \item{newdata}{data frame in which to predict --- if NULL, predictions are
for the data on which the model was fitted.  See \sQuote{Details}}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}.  See \sQuote{Details}}
  \item{type}{predictor type --- default NULL, use decomposition into
trend, signal, and noise ; other type available depending if newdata=NULL.  See \sQuote{Details}} %  "X" "TC" or "BP" in-sample predictors of Thomas-Agnan -- see reference, only if newdata=NULL
  \item{all.data}{(only applies if type="TC" and newdata is not NULL) default FALSE: return prediction only for newdata units, if TRUE return prediction for all data units.  See \sQuote{Details}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing the function to 
terminate with an error}
  \item{legacy}{(Only applies to lag and Durbin (mixed) models) default TRUE: use ad-hoc predictor, if FALSE use DGP-based predictor}
  \item{legacy.mixed}{(Only applies to mixed models) default FALSE: compute lagged variables with \eqn{[W X]_O}{[WX]o} where W =  \code{cbind(X_S, X_O)}, if TRUE compute lagged variables with \eqn{W_{OO} X_O}{Woo Xo}}
  \item{power}{(Only applies to lag and Durbin (mixed) models) use \code{powerWeights}, if default NULL, set FALSE if \code{object$method} is \dQuote{eigen}, otherwise TRUE}
  \item{order}{power series maximum limit if \code{power} is TRUE}
  \item{tol}{tolerance for convergence of power series if \code{power} is TRUE}
  \item{spChk}{should the rownames of data frames be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{x}{the object to be printed}
%   \item{pred.se}{logical, default FALSE, not available for error models}
%   \item{lagImpact}{object created by \code{impacts} with R sample draws from the fitted model}
  \item{...}{further arguments passed through}
}
\details{
In the following, the trend is the non-spatial smooth, the signal is the
spatial smooth, and the noise is the residual. The fit returned is the
sum of the trend and the signal.

The function approaches prediction first by dividing invocations between 
those with or without newdata. When no newdata is present, the response 
variable may be reconstructed as the sum of the trend, the signal, and the
noise (residuals). Since the values of the response variable are known,
their spatial lags are used to calculate signal components (Cressie 1993, p. 564). For the error
model, trend = \eqn{X \beta}{X beta}, and signal = \eqn{\lambda W y - 
\lambda W X \beta}{lambda W y - lambda W X beta}. For the lag and mixed
models, trend = \eqn{X \beta}{X beta}, and signal = \eqn{\rho W y}{rho W y}.

This approach differs from the design choices made in other software, for
example GeoDa, which does not use observations of the response variable,
and corresponds to the newdata situation described below.

When however newdata is used for prediction, no observations of the response 
variable being predicted are available. Consequently, while the trend
components are the same, the signal cannot take full account of the spatial
smooth. In the error model and Durbin error model, the signal is set to zero, since the spatial smooth is expressed in terms of the error: 
\eqn{(I - \lambda W)^{-1} \varepsilon}{inv(I - lambda W) e}.

In the lag model, the signal can be expressed in the following way (for legacy=TRUE):

\deqn{(I - \rho W) y = X \beta + \varepsilon}{(I - rho W) y = X beta + e},
\deqn{y = (I - \rho W)^{-1} X \beta + (I - \rho W)^{-1} \varepsilon}{y = inv(I - rho W) X beta + inv(I - rho W) e}

giving a feasible signal component of:

\deqn{\rho W y = \rho W (I - \rho W)^{-1} X \beta}{rho W y = rho W inv(I - rho W) X beta}

For legacy=FALSE, the trend is computed first as:

\deqn{X \beta}{X beta}

next the prediction using the DGP:

\deqn{(I - \rho W)^{-1} X \beta}{inv(I - rho W) X beta}

and the signal is found as the difference between prediction and trend. The numerical results for the legacy and DGP methods are identical.

setting the error term to zero. This also means that predictions of the
signal component for lag and mixed models require the inversion of an 
n-by-n matrix.

Because the outcomes of the spatial smooth on the error term are
unobservable, this means that the signal values for newdata are
incomplete. In the mixed model, the spatially lagged RHS variables
influence both the trend and the signal, so that the root mean square
prediction error in the examples below for this case with newdata is
smallest, although the model was not the best fit.

The all.data option is useful for some predictors, because some out-of-sample predictors can returne different prediction for in-sample units, than the same predictor type for in-sample data. 

Values for \code{type=} include \dQuote{trend}, \dQuote{TS1}, \dQuote{TC}, \dQuote{TC1}, \dQuote{BP}, \dQuote{BPW}, \dQuote{BPN}, \dQuote{KP1}, \dQuote{KP2}, \dQuote{KP3}, \dQuote{KP4}, following the notation in Thomas-Agnan et al. (2015).
}
\value{
 \code{predict.sarlm()} returns a vector of predictions with two attribute 
vectors of trend and signal values with class \code{sarlm.pred}. 
\code{print.sarlm.pred} is a print function for this class, printing and
returning a data frame with columns: "fit", "trend" and "signal".
}

\references{Haining, R. 1990 \emph{Spatial data analysis in the social and environmental sciences}, Cambridge: Cambridge University Press, p. 258; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York; Thomas-Agnan, C., Laurent, T. and Goulard, M. 2015 \emph{About predictions in spatial autoregressive models: Optimal and almost optimal strategies}, TSE Working Paper, n. 13-452.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no} and Martin Gubri}

\seealso{\code{\link{errorsarlm}}, \code{\link{lagsarlm}}}

\examples{
data(oldcol)
lw <- nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)

COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
  type="mixed")
print(p1 <- predict(COL.mix.eig))
print(p2 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, type = "trend"))
AIC(COL.mix.eig)
sqrt(deviance(COL.mix.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p1))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p2))^2)/length(COL.nb))

COL.err.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
AIC(COL.err.eig)
sqrt(deviance(COL.err.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig, newdata=COL.OLD,
  listw=lw, type = "trend")))^2)/length(COL.nb))

COL.SDerr.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
 etype="emixed")
AIC(COL.SDerr.eig)
sqrt(deviance(COL.SDerr.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig, newdata=COL.OLD,
  listw=lw, type = "trend")))^2)/length(COL.nb))

AIC(COL.lag.eig)
sqrt(deviance(COL.lag.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig, newdata=COL.OLD,
  listw=lw, type = "trend")))^2)/length(COL.nb))

p3 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, type = "trend", legacy=FALSE)
all.equal(p2, p3, check.attributes=FALSE)
p4 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, type = "trend", legacy=FALSE, power=TRUE)
all.equal(p2, p4, check.attributes=FALSE)
p5 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, type = "trend", legacy=TRUE, power=TRUE)
all.equal(p2, p5, check.attributes=FALSE)
}

\keyword{spatial}
