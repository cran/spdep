\name{asMatrixCsrListw}
\alias{asMatrixCsrListw}
\alias{asListwMatrixCsr}
\alias{asMatrixCsrIrW}
\alias{asMatrixCsrI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Convert spatial weights list to CSR matrix form}
\description{
  \code{asMatrixCsrListw} creates a \code{matrix.csr} from a \code{listw} object, \code{asMatrixCsrIrW} is a helper function used in fitting SAR models and elsewhere to make a (I - rho W) \code{matrix.csr} object from W; \code{asListwMatrixCsr} performs the reverse conversion from \code{matrix.csr} to \code{listw} object.
}
\usage{
asMatrixCsrListw(listw, zero.policy=FALSE)
asMatrixCsrIrW(W, rho)
asMatrixCsrI(n)
asListwMatrixCsr(mcsr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assume no cases of no neighbours - if they are in fact present, \code{asMatrixCsrListw()} will terminate with an error}
  \item{W}{a \code{matrix.csr} object created using \code{asMatrixCsrListw} from a \code{listw} object}
  \item{rho}{spatial regression coefficient}
  \item{n}{length of diagonal for identity matrix}
  \item{mcsr}{a \code{matrix.csr} object}
}

\value{
  \code{asMatrixCsrListw} and \code{asMatrixCsrIrW} return \code{matrix.csr} objects; \code{asListwMatrixCsr} returns a \code{listw} object
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link[SparseM]{as.matrix.csr}}}
\examples{
data(oldcol)
COL.W <- nb2listw(COL.nb, style="W")
COL.W
image(asMatrixCsrListw(COL.W))
}
\keyword{spatial}
