\name{as_dgRMatrix_listw}
\alias{as_dgRMatrix_listw}
\alias{as_dgRMatrix_listw}
\alias{as_dsTMatrix_listw}
\alias{as_dgCMatrix_I}
\alias{as_dgCMatrix_IrW}
\alias{Jacobian_W}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Interface between Matrix class objects and weights lists}
\description{
  Interface between Matrix class objects and weights lists
}
\usage{
as_dgRMatrix_listw(listw)
as_dsTMatrix_listw(listw)
as_dgCMatrix_I(n)
as_dgCMatrix_IrW(W, rho)
Jacobian_W(W, rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{W}{a \code{dsTMatrix} object created using \code{as_dsTMatrix_listw} from a symmetric \code{listw} object}
  \item{rho}{spatial regression coefficient}
  \item{n}{length of diagonal for identity matrix}
}

\value{
  Matrix package class objects
}

\author{Roger Bivand}

\examples{
example(NY_data)
\dontrun{W <- as(as_dsTMatrix_listw(listw_NY), "dsCMatrix")
I <- as(as_dgCMatrix_I(dim(W)[1]), "dtCMatrix")
IrW <- (I - 0.1 * W)
CHOL <- chol(IrW)
sum(2*log(diag(CHOL)))
sum(log(1 - 0.1 * eigenw(listw_NY)))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

