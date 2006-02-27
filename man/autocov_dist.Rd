\name{autocov_dist}
\alias{autocov_dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distance-weighted autocovariate}
\description{
  Calculates the autocovariate to be used in autonormal, autopoisson or autologistic regression. Three distance-weighting schemes are available.
}
\usage{
autocov_dist(z, xy, nbs = 1, type = "inverse", zero.policy = FALSE, style = "W", lonlat=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{the response variable}
  \item{xy}{the matrix of coordinates}
  \item{nbs}{neighbourhood radius; default is 1}
  \item{type}{the weighting scheme: "one" gives equal weight to all data points in the neighbourhood; "inverse" (the default) weights by inverse distance; "inverse.squared" weights by the square of "inverse"}
  \item{zero.policy}{If FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
  \item{style}{style' can take values W, B, C, U, and S; W gives mean values for neighbours}
  \item{lonlat}{TRUE if point coordinates are longitude-latitude decimal, in which case distances are measured in kilometers}
}

\value{
  A numeric vector of autocovariate values
}
\references{Augustin N.H., Mugglestone M.A. and Buckland S.T. (1996) An autologistic model for the spatial distribution of wildlife. \emph{Journal of Applied Ecology}, 33, 339-347; Gumpertz M.L., Graham J.M. and Ristaino J.B. (1997) Autologistic model of spatial pattern of Phytophthora epidemic in bell pepper: effects of soil variables on disease presence. \emph{Journal of Agricultural, Biological and Environmental Statistics}, 2, 131-156.}
\author{Carsten F. Dormann and Roger Bivand}

\seealso{\code{\link{nb2listw}}}
\examples{
data(columbus)
xy <- as.matrix(columbus[,13:14])
ac1a <- autocov_dist(columbus$CRIME, xy, nbs=10, style="W", type="one")
acinva <- autocov_dist(columbus$CRIME, xy, nbs=10, style="W", type="inverse")
acinv2a <- autocov_dist(columbus$CRIME, xy, nbs=10, style="W", type="inverse.squared")

plot(ac1a ~ columbus$CRIME, pch=16, asp=1)
points(acinva ~ columbus$CRIME, pch=16, col="red")
points(acinv2a ~ columbus$CRIME, pch=16, col="blue")
abline(0,1)

nb <- dnearneigh(xy, 0, 10)
lw <- nb2listw(nb, style="W")
ac1b <- lag(lw, columbus$CRIME)
all.equal(ac1b, ac1a)

nbd <- nbdists(nb, xy)
gl <- lapply(nbd, function(x) 1/x)
lw <- nb2listw(nb, glist=gl)
acinvb <- lag(lw, columbus$CRIME)
all.equal(acinvb, acinva)

gl2 <- lapply(nbd, function(x) 1/(x^2))
lw <- nb2listw(nb, glist=gl2)
acinv2b <- lag(lw, columbus$CRIME)
all.equal(acinv2b, acinv2a)

glm(CRIME ~ HOVAL + ac1b, data=columbus, family="gaussian")
spautolm(columbus$CRIME ~ HOVAL, data=columbus, listw=nb2listw(nb, style="W"))

}
\keyword{spatial}

