% Copyright 2002 by Roger S. Bivand
\name{errorsarlm}
\alias{errorsarlm}
\alias{sar.error.f}
\alias{sar.error.f.s}
\title{Spatial simultaneous autoregressive error model estimation}
\description{
  Maximum likelihood estimation of spatial simultaneous autoregressive
error models of the form:

\deqn{y = X \beta + u, u = \lambda W u + \varepsilon}{y = X beta + u, u = lambda W u + e}

where $\lambda$ is found by \code{optimize()} first and $\beta$ and other
parameters by generalized least squares subsequently.
}
\usage{
errorsarlm(formula, data=list(), listw, method="eigen", quiet=TRUE,
  zero.policy=FALSE, tol.solve=1.0e-7, tol.opt=.Machine$double.eps^0.5,
  sparsedebug=FALSE)
sar.error.f(lambda, eig, y, wy, x, WX, n, quiet)
sar.error.f.s(lambda, sn, y, wy, x, WX, n, quiet, sparsedebug)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{method}{"eigen" (default) - the Jacobian is computed as the product 
of (1 - rho*eigenvalue) using \code{eigenw}
, and "sparse" - computes the 
determinant of the sparse matrix (I - rho*W) directly using \code{log.spwdet}.
}
  \item{quiet}{default=TRUE; if FALSE, reports function values during optimization.}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{errorsarlm()} to terminate with an error}
  \item{tol.solve}{the tolerance for detecting linear dependencies in the columns of matrices to be inverted - passed to \code{solve()} (default=1.0e-7)}
  \item{tol.opt}{the desired accuracy of the optimization - passed to \code{optimize()} (default=square root of double precision machine tolerance)}
  \item{sparsedebug}{if TRUE, writes a log file on sparse matrix operations
(name sparsestats) in the current directory. To be used if sparse estimation
fails!}
  \item{lambda}{value of the spatial parameter}
  \item{eig}{eigenvalues of the full spatial weights matrix from \code{eigenw}}
  \item{y}{dependent variable}
  \item{wy}{spatially lagged dependent variable}
  \item{x}{independent variables}
  \item{WX}{spatially lagged independent variables}
  \item{n}{length of y (and eig)}
  \item{sn}{sparse spatial neighbour object from \code{listw2sn}}
}
\details{
  The asymptotic standard error of $\lambda$ is only computed when 
method=eigen, because the full matrix operations involved would be costly
for large n typically associated with the choice of method=sparse. 

When using the sparse method, the user takes (unfortunately) full
responsibility for possible failures, including R terminating with
a core dump! Safeguards have been put in place to try to trap errant
behaviour in the sparse functions' memory allocation, but they may
not always help. When sparsedebug is TRUE, a log file (sparsestats)
is written in the working directory - the figure of interest is the
number of allocated blocks. At present, \code{spwdet} will fail when
this increases over the number initially allocated, but will not release
memory allocated by the sparse functions. In the event of problems,
save your workspace and quit R. Problems seem to be related to larger
n, and to an unknown trigger precipitating incontrolled fillin, in the
course of which the sparse routines lose track of their memory pointers,
and then provoke a segmentation fault trying to free unallocated memory.
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

The sar.error.* functions return the value of the log likelihood function
at $\lambda$.
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
\code{\link{eigenw}}, \code{\link{log.spwdet}}, \code{\link{predict.sarlm}},
\code{\link{residuals.sarlm}}
}

\examples{
data(oldcol)
COL.err.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="eigen", quiet=FALSE)
COL.err.sp <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="sparse", quiet=FALSE)
summary(COL.err.eig)
summary(COL.err.sp)
}
\keyword{spatial}
