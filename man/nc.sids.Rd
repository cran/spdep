\name{nc.sids}
\alias{nc.sids}
\alias{nc.utmbbs}
\alias{nc.utm.polys}
\alias{sidsorig.nb}
\alias{utm18.countyseats}
\non_function{}
\title{North Carolina SIDS data}
\usage{data(nc.sids)}
\description{
The \code{nc.sids} data frame has 100 rows and 9 columns. It contains data
given in Cressie (1991, pp. 386-9) and elsewhere on sudden infant deaths in 
North Carolina for 1974-78 and 1979-84. The data set also contains polygon
boundaries for the counties reprojected to UTM zone 18, the neighbour list 
given by Cressie omitting self-neighbours, and the coordinates of the county 
seats approximately reprojected to UTM zone 18.

}
\format{
  This data frame contains the following columns:
  \describe{
    \item{CNTY.ID}{county ID}
    \item{BIR74}{births, 1974-78}
    \item{SID74}{SID deaths, 1974-78}
    \item{NWBIR74}{non-white births, 1974-78}
    \item{BIR79}{births, 1979-84}
    \item{SID79}{SID deaths, 1979-84}
    \item{NWBIR79}{non-white births, 1979-84}
    \item{east}{eastings, county seat, miles, local projection}
    \item{north}{northings, county seat, miles, local projection}
  }
}
\source{
 Cressie, N (1991), Statistics for spatial data. New York: Wiley, pp. 386-389.
 \url{http://www.spacestat.com/data4.htm}.
}

\examples{
data(nc.sids)
stem(floor(sort(1000*(nc.sids$SID74 + 1)/nc.sids$BIR74)*10)/10, scale=2)
sids.phat <- sum(nc.sids$SID74) / sum(nc.sids$BIR74)
pm <- ppois(nc.sids$SID74, sids.phat*nc.sids$BIR74)
pm.f <- as.ordered(cut(pm, breaks=c(0.0, 0.01, 0.05, 0.1, 0.9, 0.95, 0.99, 1),
  include.lowest=TRUE))
cols <- cm.colors(7)
plotpolys(nc.utm.polys, nc.utmbbs, col=cols[codes(pm.f)])
legend(c(-280, -70), c(3700, 3900), legend=paste("prob.", levels(pm.f)),
  fill=cols, bty="n")
ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) +
  sqrt((nc.sids$SID74+1)/nc.sids$BIR74))
names(ft.SID74) <- rownames(nc.sids)
stem(ft.SID74, scale=2)
plotpolys(nc.utm.polys, nc.utmbbs, border="grey")
plot(sidsorig.nb, utm18.countyseats, add=TRUE, col="blue")
text(utm18.countyseats[card(sidsorig.nb) == 0,],
  rownames(nc.sids)[card(sidsorig.nb) == 0], pos=3)
moran.plot(ft.SID74, nb2listw(sidsorig.nb, zero.policy=TRUE),
  zero.policy=TRUE)
drop.no.neighs <- !(1:length(sidsorig.nb) \%in\% which(card(sidsorig.nb) == 0))
sub.sidsorig.nb <- subset(sidsorig.nb, drop.no.neighs)
sub.ft.SID74 <- subset(ft.SID74, drop.no.neighs)
moran.plot(sub.ft.SID74, nb2listw(sub.sidsorig.nb),
  labels=attr(sub.sidsorig.nb, "region.id"))
tr.SID74 <- ft.SID74*sqrt(nc.sids$BIR74)
names(tr.SID74) <- names(ft.SID74)
moran.test(tr.SID74, nb2listw(sidsorig.nb, zero.policy=TRUE),
 zero.policy=TRUE, alternative="two.sided")
geary.test(tr.SID74, nb2listw(sidsorig.nb, zero.policy=TRUE),
 zero.policy=TRUE, alternative="two.sided")
}
\keyword{datasets}
