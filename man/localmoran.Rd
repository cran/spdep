% Copyright 2001 by Roger S. Bivand
\name{localmoran}
\alias{localmoran}
\title{Local Moran's I statistic}
\description{
 The local spatial statistic Moran's I is calculated for each zone based on the
spatial weights object used. The values returned include a Z-value, and may be
used as a diagnostic tool. The statistic is:
\deqn{I_i = \frac{(x_i-\bar{x})}{{\sum_{k=1}^{n}(x_k-\bar{x})^2}/n}{\sum_{j=1}^{n}w_{ij}(x_j-\bar{x})}},
and its expectation and variance are given in Anselin (1995).
}
\usage{
localmoran(x, listw, zero.policy=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
}
\value{

  \item{Ii}{local moran statistic}
  \item{E.Ii}{expectation of local moran statistic}
  \item{Var.Ii}{variance of local moran statistic}
  \item{Z.Ii}{standard deviate of local moran statistic}
}
\references{Anselin, L. 1995. Local indicators of spatial association,
Geographical Analysis, 27, 93--115;
Getis, A. and Ord, J. K. 1996 Local spatial
statistics: an overview. In P. Longley and M. Batty (eds) \emph{Spatial
analysis: modelling in a GIS environment} (Cambridge: Geoinformation
International), 261--277.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{localG}}}

\examples{
data(afcon)
oid <- order(afcon$id)
resG <- localG(afcon$totcon, nb2listw(include.self(paper.nb)))
print(data.frame(resG[oid], row.names=afcon$name[oid]), digits=2)
resI <- localmoran(afcon$totcon, nb2listw(paper.nb))
print(data.frame(resI[oid,], row.names=afcon$name[oid]), digits=2)

}
\keyword{spatial}
