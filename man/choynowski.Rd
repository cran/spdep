% Copyright 2004 by Roger S. Bivand
\name{choynowski}
\alias{choynowski}
\title{Choynowski probability map values}
\description{
  Calculates Choynowski probability map values.
}
\usage{
choynowski(n, x, row.names=NULL, tol = .Machine$double.eps^0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{a numeric vector of counts of cases}
  \item{x}{a numeric vector of populations at risk}
  \item{row.names}{row names passed through to output data frame}
  \item{tol}{accumulate values for observed counts >= expected until value less than tol}
}

\value{
  A data frame with columns:
  \item{pmap}{Poisson probability map values: probablility of getting a more ``extreme'' count than actually observed, one-tailed with less than expected and more than expected folded together}
  \item{type}{logical: TRUE if observed count less than expected}
}
\references{Choynowski, M (1959) Maps based on probabilities, Journal of
the American Statistical Association, 54, 385--388; Cressie, N, Read,
TRC (1985), Do sudden infant deaths come in clusters? Statistics and
Decisions, Supplement Issue 2, 333--349; Bailey T, Gatrell A (1995)
Interactive Spatial Data Analysis, Harlow: Longman, pp. 300--303.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{probmap}}}
\examples{
data(auckland)
res <- choynowski(auckland$Deaths.1977.85, 9*auckland$Under.5.1981)
res1 <- probmap(auckland$Deaths.1977.85, 9*auckland$Under.5.1981)
table(abs(res$pmap - res1$pmap) < 0.00001, res$type)
plot(auckpolys, forcefill=FALSE)
lt005 <- (res$pmap < 0.05) & (res$type)
ge005 <- (res$pmap < 0.05) & (!res$type)
plot(subset(auckpolys, lt005), add=TRUE, col=grey(2/7), forcefill=FALSE) 
plot(subset(auckpolys, ge005), add=TRUE, col=grey(5/7), forcefill=FALSE)
legend(c(70,90), c(70,95), fill=grey(c(2,5)/7),
 legend=c("low", "high"), bty="n")
}
\keyword{spatial}
