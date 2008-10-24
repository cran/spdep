\name{aple.mc}
\alias{aple.mc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Approximate profile-likelihood estimator (APLE) permutation test}
\description{
  A permutation bootstrap test for the approximate profile-likelihood estimator (APLE).
}
\usage{
aple.mc(x, listw, nsim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a zero-mean detrended continuous variable}
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{nsim}{number of simulations}
}

\value{
  A \code{boot} object as reurned by the \code{boot} function.
}
\references{Li, H, Calder, C. A. and Cressie N. A. C. (2007) Beyond Moran's I: testing for spatial dependence based on the spatial autoregressive model. Geographical Analysis 39, 357-375.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{aple}}, \code{\link[boot]{boot}}}
\examples{
example(aple)
boot_out <- aple.mc(scale(wheat$yield_detrend, scale=FALSE),
 nb2listw(nbr12, style="W"), nsim=500)
plot(boot_out)
boot_out
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

