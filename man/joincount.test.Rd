% Copyright 2001 by Roger S. Bivand
\name{joincount.test}
\alias{joincount.test}
\alias{print.jclist}
\title{BB join count statistic for k-coloured factors}
\description{
The BB join count test for spatial autocorrelation using a spatial
weights matrix in weights list form for testing whether same-colour joins
occur more frequently than would be expected if the zones were labelled
in a spatially random way. The assumptions underlying the test are
sensitive to the form of the graph of neighbour relationships and other
factors, and results may be checked against those of \code{geary.mc}
permutations.
}
\usage{
joincount.test(fx, listw, alternative="greater")
print.jclist(x, ...)
}
\arguments{
  \item{fx}{a factor of the same length as the neighbours and weights objects in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided.}
  \item{x}{object to be printed}
  \item{...}{arguments to be passed through for printing}
}

}
\value{
A list with class \code{jclist} of lists with class \code{htest} for each of the k colours containing the following components:
  \item{statistic}{the value of the standard deviate of the join count statistic.}
  \item{p.value}{the p-value of the test.}
  \item{estimate}{the value of the observed statistic, its expectation and variance.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data.}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 20.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{Probability values one-tailed, as the hypothesis most often used is no autocorrelation against positive autocorrelation.}

\seealso{\code{\link{joincount}}, \code{\link{joincount.mc}}}

\examples{
data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
joincount.test(HICRIME, nb2listw(COL.nb, style="B"))
joincount.test(HICRIME, nb2listw(COL.nb, style="C"))
joincount.test(HICRIME, nb2listw(COL.nb, style="S"))
joincount.test(HICRIME, nb2listw(COL.nb, style="W"))
by(card(COL.nb), HICRIME, summary)
}
\keyword{spatial}
