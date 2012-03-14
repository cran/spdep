% Copyright 2002 by Roger S. Bivand
\name{predict.sarlm}
\alias{predict.sarlm}
\alias{print.sarlm.pred}

\title{Prediction for spatial simultaneous autoregressive linear
model objects}
\description{
  \code{predict.sarlm()} calculates predictions as far as is at present
possible for for spatial simultaneous autoregressive linear
model objects, using Haining's terminology for decomposition into
trend, signal, and noise --- see reference.
}
\usage{
\method{predict}{sarlm}(object, newdata = NULL, listw = NULL,
 zero.policy = NULL, legacy=TRUE, power=NULL, order=250,
 tol=.Machine$double.eps^(3/5), pred.se=FALSE, lagImpact=NULL, ...)
\method{print}{sarlm.pred}(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{object}{\code{sarlm} object returned by \code{lagsarlm} or 
\code{errorsarlm}}
  \item{newdata}{Data frame in which to predict --- if NULL, predictions are
for the data on which the model was fitted}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{default NULL, use global option value; if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing the function to 
terminate with an error}
  \item{legacy}{(Only applies to lag and Durbin (mixed) models) default TRUE: use ad-hoc predictor, if FALSE use DGP-based predictor}
  \item{power}{(Only applies to lag and Durbin (mixed) models) use \code{powerWeights}, if default NULL, set FALSE if \code{object$method} is \dQuote{eigen}, otherwise TRUE}
  \item{order}{Power series maximum limit if \code{power} is TRUE}
   \item{tol}{Tolerance for convergence of power series if \code{power} is TRUE}
   \item{x}{the object to be printed}
   \item{pred.se}{logical, default FALSE, not available for error models}
   \item{lagImpact}{object created by \code{impacts} with R sample draws from the fitted model}
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
smallest, although the model was not the best fit 
}
\value{
 \code{predict.sarlm()} returns a vector of predictions with two attribute 
vectors of trend and signal values with class \code{sarlm.pred}. 
\code{print.sarlm.pred} is a print function for this class, printing and
returning a data frame with columns: "fit", "trend" and "signal".
}

\references{Haining, R. 1990 \emph{Spatial data analysis in the social
and environmental sciences}, Cambridge: Cambridge University Press, p. 258; Cressie, N. A. C. 1993 \emph{Statistics for spatial data}, Wiley, New York.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{errorsarlm}}, \code{\link{lagsarlm}}}

\examples{
data(oldcol)
lw <- nb2listw(COL.nb)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
  type="mixed")
COL.err.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw)
COL.SDerr.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, lw,
 etype="emixed")
print(p1 <- predict(COL.mix.eig))
print(p2 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw))
AIC(COL.mix.eig)
sqrt(deviance(COL.mix.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p1))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p2))^2)/length(COL.nb))
AIC(COL.err.eig)
sqrt(deviance(COL.err.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig, newdata=COL.OLD,
  listw=lw)))^2)/length(COL.nb))
AIC(COL.SDerr.eig)
sqrt(deviance(COL.SDerr.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.SDerr.eig, newdata=COL.OLD,
  listw=lw)))^2)/length(COL.nb))
AIC(COL.lag.eig)
sqrt(deviance(COL.lag.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig, newdata=COL.OLD,
  listw=lw)))^2)/length(COL.nb))
p3 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, legacy=FALSE)
all.equal(p2, p3)
p4 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, legacy=FALSE, power=TRUE)
all.equal(p2, p4)
p5 <- predict(COL.mix.eig, newdata=COL.OLD, listw=lw, legacy=TRUE, power=TRUE)
all.equal(p2, p5)
}

\keyword{spatial}
