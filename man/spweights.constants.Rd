% Copyright 2001 by Roger S. Bivand
\name{spweights.constants}
\alias{spweights.constants}
\alias{Szero}
\title{Provides constants for spatial weights matrices}
\description{
 The function calculates the constants needed for tests of spatial autocorrelation for general weights matrices represented as \code{listw} objects
}
\usage{
spweights.constants(listw, zero.policy=FALSE)
Szero(listw)
}
\arguments{
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{zero.policy}{if TRUE ignore zones without neighbours, if FALSE fail when encountered}
}

\value{

  \item{n}{number of zones}
  \item{n1}{n - 1}
  \item{n2}{n - 2}
  \item{n3}{n - 3}
  \item{nn}{n * n}
  \item{S0}{global sum of weights}
  \item{S1}{S1 sum of weights}
  \item{S2}{S2 sum of weights}
}
\references{Haining, R. 1990 Spatial data analysis in the social and environmental sciences, Cambridge University Press, p. 233; Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 19, 21.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link[spweights]{nb2listw}}}

\examples{
data(oldcol)
B <- spweights.constants(nb2listw(COL.nb, style="B"))
W <- spweights.constants(nb2listw(COL.nb, style="W"))
C <- spweights.constants(nb2listw(COL.nb, style="C"))
S <- spweights.constants(nb2listw(COL.nb, style="S"))
U <- spweights.constants(nb2listw(COL.nb, style="U"))
print(data.frame(rbind(unlist(B), unlist(W), unlist(C), unlist(S), unlist(U)),
  row.names=c("B", "W", "C", "S", "U")))
}
\keyword{spatial}
