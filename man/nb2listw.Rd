% Copyright 2001 by Roger S. Bivand
\name{nb2listw}
\alias{nb2listw}
\title{Spatial weights for neighbours lists}
\description{
The function supplements a neighbours list with spatial weights for the chosen coding scheme.
}
\usage{
nb2listw(neighbours, glist=NULL, style="W", zero.policy=FALSE)
}
\arguments{
  \item{neighbours}{an object of class \code{nb}}
  \item{glist}{list of general weights corresponding to neighbours}
  \item{style}{\code{style} can take values W, B, C, and S}
  \item{zero.policy}{If FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
}
\details{
Starting from a binary neighbours list, in which regions are either listed as neighbours or are absent (thus not in the set of neighbours for some definition), the function adds a weights list with values given by the coding scheme style chosen. B is the basic binary coding, W is row standardised, C is globally standardised, while S is the variance-stabilizing coding scheme proposed by Tiefelsdorf et al. 1999, p. 167-168.

If zero policy is set to TRUE, weights vectors of zero length are inserted for regions without neighbour in the neighbours list. These will in turn generate lag values of zero, equivalent to the sum of products of the zero row \code{t(rep(0, length=length(neighbours))) \%*\% x}, for arbitraty numerical vector \code{x} of length \code{length(neighbours)}. The spatially lagged value of x for the zero-neighbour region will then be zero, which may (or may not) be a sensible choice.
}
\value{
A \code{listw} object with the following members:
  \item{style}{one of W, B, C, S as above}
  \item{neighbours}{the input neighbours list}
  \item{weights}{the weights for the neighbours and chosen style, with attributes set to report the type of relationships (binary or general, if general the form of the glist argument), and style as above}
}
\references{Tiefelsdorf, M., Griffith, D. A., Boots, B. 1999 A variance-stabilizing coding scheme for spatial link matrices, Environment and Planning A, 31, pp. 165-180.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{summary.nb}}, \code{\link{read.gal}}}

\examples{
data(columbus)
cards <- card(col.gal.nb)
col.w <- nb2listw(col.gal.nb)
plot(cards, unlist(lapply(col.w$weights, sum)),xlim=c(0,10),
ylim=c(0,10), xlab="number of links", ylab="row sums of weights")
col.b <- nb2listw(col.gal.nb, style="B")
points(cards, unlist(lapply(col.b$weights, sum)), col="red")
col.c <- nb2listw(col.gal.nb, style="C")
points(cards, unlist(lapply(col.c$weights, sum)), col="green")
col.s <- nb2listw(col.gal.nb, style="S")
points(cards, unlist(lapply(col.s$weights, sum)), col="blue")
legend(x=c(0, 1), y=c(7, 9), legend=c("W", "B", "C", "S"),
col=c("black", "red", "green", "blue"), pch=rep(1,4))
dlist <- nbdists(col.gal.nb, coords)
dlist <- lapply(dlist, function(x) 1/x)
col.w.d <- nb2listw(col.gal.nb, glist=dlist)
summary(unlist(col.w$weights))
summary(unlist(col.w.d$weights))
}

\keyword{spatial}
