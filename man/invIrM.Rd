% Copyright 2001 by Roger S. Bivand
\name{invIrM}
\alias{invIrM}

\title{Compute SAR generating operator}
\description{
  Computes the matrix used for generating simultaneous autoregressive random variables, for a given value of rho, a neighbours list object, a chosen coding scheme style, and optionally a list of general weights corresponding to neighbours.
}
\usage{
invIrM(neighbours, rho, glist=NULL, style="W")
}

\arguments{
  \item{neighbours}{an object of class \code{nb}}
  \item{rho}{autoregressive parameter}
  \item{glist}{list of general weights corresponding to neighbours}
  \item{style}{\code{style} can take values W, B, C, and S}
}
\details{
The function generates the full weights matrix V, checks that rho lies in its feasible range between 1/min(eigen(V)) and 1/max(eigen(V)), and returns the nxn inverted matrix \deqn{(I - \rho V)^{-1}}
  
}
\value{
  An nxn matrix with a "call" attribute.
}
\references{Tiefelsdorf, M., Griffith, D. A., Boots, B. 1999 A variance-stabilizing coding scheme for spatial link matrices, Environment and Planning A, 31, pp. 165-180; Tiefelsdorf, M. 2000 Modelling spatial processes, Lecture notes in earth sciences, Springer, p. 76; Haining, R. 1990 Spatial data analysis in the social and environmental sciences, Cambridge University Press, p. 117; Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 152.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}}

\examples{

nb7rt <- cell2nb(7, 7, torus=TRUE)
x <- matrix(rnorm(500*length(nb7rt)), nrow=length(nb7rt))
res0 <- apply(invIrM(nb7rt, rho=0.0) \%*\% x, 2, function(x) var(x)/length(x))
res2 <- apply(invIrM(nb7rt, rho=0.2) \%*\% x, 2, function(x) var(x)/length(x))
res4 <- apply(invIrM(nb7rt, rho=0.4) \%*\% x, 2, function(x) var(x)/length(x))
res6 <- apply(invIrM(nb7rt, rho=0.6) \%*\% x, 2, function(x) var(x)/length(x))
res8 <- apply(invIrM(nb7rt, rho=0.8) \%*\% x, 2, function(x) var(x)/length(x))
res9 <- apply(invIrM(nb7rt, rho=0.9) \%*\% x, 2, function(x) var(x)/length(x))
plot(density(res9), col="red", xlim=c(-0.01, max(density(res9)$x)),
  ylim=range(density(res0)$y),
  xlab="estimated variance of the mean",
  main=expression(paste("Effects of spatial autocorrelation for different ",
    rho, " values")))
lines(density(res0), col="black")
lines(density(res2), col="brown")
lines(density(res4), col="green")
lines(density(res6), col="orange")
lines(density(res8), col="pink")
legend(c(-0.02, 0.01), c(7, 25), legend=c("0.0", "0.2", "0.4", "0.6", "0.8", "0.9"), col=c("black", "brown", "green", "orange", "pink", "red"), lty=1, bty="n")

}
\keyword{spatial}
