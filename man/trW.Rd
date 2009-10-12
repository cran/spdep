\name{trW}
\Rdversion{1.1}
\alias{trW}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial weights matrix powers traces}
\description{The function is used to prepare a vector of traces of pwers of a spatial weights matrix}
\usage{
trW(W, m = 100, p = 50, type = "mult")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{W}{A spatial weights matrix in CsparseMatrix form}
  \item{m}{The number of powers}
  \item{p}{The number of samples used in Monte Carlo simulation of the traces if type is MC}
  \item{type}{Either mult (default) for powering a sparse matrix, or MC for Monte Carlo simulation of the traces}
}

\value{
A numeric vector of \code{m} traces.
}
\references{LeSage J and RK Pace (2009) Introduction to Spatial Econometrics. CRC Press, Boca Raton, pp. 96--105.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}


\seealso{
\code{\link{as_dgRMatrix_listw}}, \code{\link{nb2listw}}
}
\examples{
example(columbus)
listw <- nb2listw(col.gal.nb)
W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
trMat <- trW(W, type="mult")
str(trMat)
trMC <- trW(W, type="MC")
str(trMC)
plot(trMat, trMC)
}
\keyword{spatial}

