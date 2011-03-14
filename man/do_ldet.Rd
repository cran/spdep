\name{do_ldet}
\alias{do_ldet}
\alias{eigen_setup}
\alias{mcdet_setup}
\alias{cheb_setup}
\alias{spam_setup}
\alias{spam_update_setup}
\alias{Matrix_setup}
\alias{Matrix_J_setup}
\alias{LU_setup}
\alias{moments_setup}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial regression model Jacobian computations}

\description{These functions are made available in the package namespace for other developers, and are not intended for users. They provide a shared infrastructure for setting up data for Jacobian computation, and then for caclulating the Jacobian, either exactly or approximately, in maximum likelihood fitting of spatial regression models. The techniques used are the exact eigenvalue, Cholesky decompositions (Matrix, spam), and LU ones, with Chebyshev and Monte Carlo approximations; moments use the methods due to Martin and Smirnov/Anselin.}

\usage{
do_ldet(coef, env, which=1)
cheb_setup(env, q=5, which=1)
mcdet_setup(env, p=16, m=30, which=1)
eigen_setup(env, which=1)
spam_setup(env, pivot="MMD", which=1)
spam_update_setup(env, in_coef=0.1, pivot="MMD", which=1)
Matrix_setup(env, Imult, super=as.logical(NA), which=1)
Matrix_J_setup(env, super=FALSE, which=1)
LU_setup(env, which=1)
moments_setup(env, trs=NULL, m, p, type="MC", correct=TRUE, trunc=TRUE, which=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{coef}{spatial coefficient value}
  \item{env}{environment containing pre-computed objects, fixed after assignment in setup functions}
  \item{which}{default 1; if 2, use second listw object}
  \item{q}{Chebyshev approximation order; default in calling spdep functions is 5, here it cannot be missing and does not have a default}
  \item{p}{Monte Carlo approximation number of random normal variables; default calling spdep functions is 16, here it cannot be missing and does not have a default}
  \item{m}{Monte Carlo approximation number of series terms; default in calling spdep functions is 30, here it cannot be missing and does not have a default; \code{m} serves the same purpose in the moments method}
  \item{pivot}{default \dQuote{MMD}, may also be \dQuote{RCM} for Cholesky decompisition using spam}
  \item{in_coef}{fill-in initiation coefficient value, default 0.1}
  \item{Imult}{see \code{\link[Matrix]{Cholesky}}; numeric scalar which defaults to zero. The matrix that is decomposed is A+m*I where m is the value of Imult and I is the identity matrix of order ncol(A). Default in calling spdep functions is 2, here it cannot be missing and does not have a default, but is rescaled for binary weights matrices in proportion to the maximim row sum in those calling functions}
  \item{super}{see \code{\link[Matrix]{Cholesky}}; logical scalar indicating is a supernodal decomposition should be created.  The alternative is a simplicial decomposition. Default in calling spdep functions is FALSE for \dQuote{Matrix_J} and \code{as.logical(NA)} for \dQuote{Matrix}.  Setting it to NA leaves the choice to a CHOLMOD-internal heuristic}
  \item{trs}{A numeric vector of \code{m} traces, as from \code{trW}}
  \item{type}{moments trace type, see \code{\link{trW}}}
  \item{correct}{default TRUE: use Smirnov correction term, see \code{\link{trW}}}
  \item{trunc}{default TRUE: truncate Smirnov correction term, see \code{\link{trW}}}
}

\details{Since environments are containers in the R workspace passed by reference rather than by value, they are useful for passing objects to functions called in numerical optimisation, here for the maximum likelihood estimation of spatial regression models. This technique can save a little time on each function call, balanced against the need to access the objects in the environment inside the function. The environment should contain a \code{family} string object either \dQuote{SAR}, \dQuote{CAR} or \dQuote{SMA} (used in \code{do_ldet} to choose spatial moving average in \code{spautolm}, and these specific objects before calling the set-up functions:

\describe{
  \item{eigen}{Classical Ord eigenvalue computations:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
    \item{verbose}{logical scalar: legacy report print control, for historical reasons only}
  } and assigns to the environment:
  \describe{
    \item{eig}{a vector of eigenvalues}
    \item{eig.range}{the search interval for the spatial coefficient}
    \item{method}{string: \dQuote{eigen}}
  }
}
  \item{Matrix}{Sparse matrix pre-computed Cholesky decomposition with fast updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
  } and assigns to the environment:
  \describe{
    \item{csrw}{sparse spatial weights matrix}
    \item{nW}{negative sparse spatial weights matrix}
    \item{pChol}{a \dQuote{CHMfactor} from factorising \code{csrw} with \code{\link[Matrix]{Cholesky}}}
    \item{nChol}{a \dQuote{CHMfactor} from factorising \code{nW} with \code{\link[Matrix]{Cholesky}}}
    \item{method}{string: \dQuote{Matrix}}
  }
}
  \item{Matrix_J}{Standard Cholesky decomposition without updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{csrw}{sparse spatial weights matrix}
    \item{I}{sparse identity matrix}
    \item{super}{the value of the \code{super} argument}
    \item{method}{string: \dQuote{Matrix_J}}
  }
}
  \item{spam}{Standard Cholesky decomposition without updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{csrw}{sparse spatial weights matrix}
    \item{I}{sparse identity matrix}
    \item{pivot}{string --- pivot method}
    \item{method}{string: \dQuote{spam}}
  }
}
  \item{spam_update}{Pre-computed Cholesky decomposition with updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{csrw}{sparse spatial weights matrix}
    \item{I}{sparse identity matrix}
    \item{csrwchol}{A Cholesky decomposition for updating}
    \item{method}{string: \dQuote{spam}}
  }
}
  \item{LU}{Standard LU decomposition without updating:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{n}{number of spatial objects}
  } and assigns to the environment:
  \describe{
    \item{W}{sparse spatial weights matrix}
    \item{I}{sparse identity matrix}
    \item{method}{string: \dQuote{LU}}
  }
}
  \item{MC}{Monte Carlo approximation:
  \describe{
    \item{listw}{A listw spatial weights object}
  } and assigns to the environment:
  \describe{
    \item{clx}{list of Monte Carlo approximation terms}
    \item{W}{sparse spatial weights matrix}
    \item{method}{string: \dQuote{MC}}
  }
}
  \item{cheb}{Chebyshev approximation:
  \describe{
    \item{listw}{A listw spatial weights object}
  } and assigns to the environment:
  \describe{
    \item{trT}{vector of Chebyshev approximation terms}
    \item{W}{sparse spatial weights matrix}
    \item{method}{string: \dQuote{Chebyshev}}
  }
}
  \item{moments}{moments approximation:
  \describe{
    \item{listw}{A listw spatial weights object}
    \item{can.sim}{logical scalar: can the spatial weights be made symmetric by similarity}
  } and assigns to the environment:
  \describe{
    \item{trs}{vector of Chebyshev approximation terms}
    \item{correct}{logical scalar: use Smirnov correction term}
    \item{trunc}{logical scalar: truncate Smirnov correction term}
    \item{method}{string: \dQuote{moments}}
  }
}
}

Some set-up functions may also assign \code{similar} to the environment if the weights were made symmetric by similarity.
}

\value{\code{do_ldet} returns the value of the Jacobian for the calculation method recorded in the environment argument; the remaining functions modify the environment in place as a side effect and return nothing.}

\references{LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton, pp. 77--110}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{spautolm}}, \code{\link{lagsarlm}}, \code{\link{errorsarlm}}, \code{\link[Matrix]{Cholesky}}}

\examples{
data(boston)
lw <- nb2listw(boston.soi)
can.sim <- spdep:::can.be.simmed(lw)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("verbose", FALSE, envir=env)
assign("family", "SAR", envir=env)
eigen_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
assign("n", length(boston.soi), envir=env)
Matrix_setup(env, Imult=2, super=FALSE)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
spam_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
LU_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
cheb_setup(env, q=5)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
set.seed(12345)
mcdet_setup(env, p=16, m=30)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

