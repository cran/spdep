% Copyright 2003 by Roger S. Bivand
\name{probmap}
\alias{probmap}
\title{Probability mapping for rates}
\description{
  The function returns a data frame of rates for counts in populations at risk with crude rates, expected counts of cases, relative risks, and Poisson probabilities.
}
\usage{
probmap(n, x)
}
\arguments{
  \item{n}{a numeric vector of counts of cases}
  \item{x}{a numeric vector of populations at risk}
}
\details{
  The function returns a data frame, from which rates may be mapped after class intervals have been chosen. The class intervals used in the examples are mostly taken from the referenced source.
}
\value{
  \item{raw}{raw (crude) rates}
  \item{expCount}{expected counts of cases assuming global rate}
  \item{relRisk}{relative risks: ratio of observed and expected counts of cases multiplied by 100}
  \item{pmap}{Poisson probability map values: probablility of getting a more ``extreme'' count than actually observed - here two-tailed, with extreme tails indicating ``unusual'' values}
}
\references{Bailey T, Gatrell A (1995) Interactive Spatial Data Analysis, Harlow: Longman, pp. 300--303.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{EBest}}, \code{\link{EBlocal}}, \code{\link{ppois}}}

\examples{
data(auckland)
res <- probmap(auckland$Deaths.1977.85, 9*auckland$Under.5.1981)
res.f <- as.ordered(cut(res$raw*1000, breaks=c(-Inf,2,2.5,3,3.5,Inf),
 include.lowest=TRUE))
cols <- grey(6:2/7)
plotpolys(auckpolys, auckbbs, col=cols[codes(res.f)])
legend(c(70,90), c(70,95), fill=cols, legend=levels(res.f), bty="n")
title(main="Crude (raw) estimates of infant mortality per 1000 per year")
res.f <- as.ordered(cut(res$relRisk, breaks=c(-Inf,47,83,118,154,190,Inf),
 include.lowest=TRUE))
cols <- cm.colors(6)
plotpolys(auckpolys, auckbbs, col=cols[codes(res.f)])
legend(c(70,90), c(70,95), fill=cols, legend=levels(res.f), bty="n")
title(main="Standardised mortality ratios for Auckland child deaths")
res.f <- as.ordered(cut(res$pmap, breaks=c(0,0.05,0.1,0.2,0.8,0.9,0.95,1),
 include.lowest=TRUE))
cols <- cm.colors(7)
plotpolys(auckpolys, auckbbs, col=cols[codes(res.f)])
legend(c(70,90), c(70,95), fill=cols, legend=levels(res.f), bty="n")
title(main="Poisson probabilities for Auckland child mortality")
}
\keyword{spatial}