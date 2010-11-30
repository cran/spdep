\name{trW}
\Rdversion{1.1}
\alias{trW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial weights matrix powers traces}
\description{The function is used to prepare a vector of traces of pwers of a spatial weights matrix}
\usage{
trW(W=NULL, m = 30, p = 16, type = "mult", listw=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{A spatial weights matrix in CsparseMatrix form}
  \item{m}{The number of powers; must be an even number for \sQuote{type}=\dQuote{moments} (default changed from 100 to 30 (2010-11-17))}
  \item{p}{The number of samples used in Monte Carlo simulation of the traces if type is MC (default changed from 50 to 16 (2010-11-17))}
  \item{type}{Either \dQuote{mult} (default) for powering a sparse matrix (with moderate or larger N, the matrix becomes dense, and may lead to swapping), or \dQuote{MC} for Monte Carlo simulation of the traces, or \dQuote{moments} to use the looping space saving algorithm proposed by Smirnov and Anselin (2009) - for \dQuote{moments}, \code{W} must be symmetric, for row-standardised weights through a similarity transformation}
  \item{listw}{a listw object, which should either be fully symmetric, or be constructed as similar to symmetric from intrinsically symmetric neighbours using \code{\link{similar.listw}}, used with \sQuote{type}=\dQuote{moments}}
}

\value{
A numeric vector of \code{m} traces.
}
\references{LeSage J and RK Pace (2009) \emph{Introduction to Spatial Econometrics}. CRC Press, Boca Raton, pp. 96--105; Smirnov O and L Anselin (2009) An O(N) parallel method of computing the Log-Jacobian of the variable transformation for models with spatial interaction on a lattice. \emph{Computational Statistics and Data Analysis} 53 (2009) 2983--2984.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{
\code{\link{as_dgRMatrix_listw}}, \code{\link{nb2listw}}
}
\examples{
example(columbus)
listw <- nb2listw(col.gal.nb)
W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
system.time(trMat <- trW(W, type="mult"))
str(trMat)
set.seed(1)
system.time(trMC <- trW(W, type="MC"))
str(trMC)
plot(trMat, trMC)
listwS <- similar.listw(listw)
W <- as(as(as_dgRMatrix_listw(listwS), "CsparseMatrix"), "symmetricMatrix")
system.time(trmom <- trW(W, m=24, type="moments"))
str(trmom)
all.equal(trMat[1:24], trmom, check.attributes=FALSE)
system.time(trMat <- trW(W, m=24, type="mult"))
str(trMat)
all.equal(trMat, trmom, check.attributes=FALSE)
set.seed(1)
system.time(trMC <- trW(W, m=24, type="MC"))
str(trMC)
}
\keyword{spatial}

