\name{nc.sids}
\alias{nc.sids}
\alias{ncCR85.nb}
\alias{ncCC89.nb}
\non_function{}
\title{North Carolina SIDS data}
\usage{data(nc.sids)}
\description{
The \code{nc.sids} data frame has 100 rows and 9 columns. It contains data
given in Cressie (1991, pp. 386-9), Cressie and Read (1985) and Cressie and Chan (1989) on sudden infant deaths in North Carolina for 1974-78 and 1979-84. The data set also contains the neighbour list given by Cressie and Chan (1989) omitting self-neighbours (ncCC89.nb), and the neighbour list given by Cressie and Read (1985) for contiguities (ncCR85.nb). The data are ordered by county ID number, not alphabetically as in the source tables.

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
    \item{x}{easting, county seats approximately reprojected to UTM zone 18}
    \item{y}{northings, county seats approximately reprojected to UTM zone 18}
    \item{lon}{easting, county seats in long-lat}
    \item{lat}{northings, county seats in long-lat}
    \item{L.id}{Cressie and Read (1985) L index}
    \item{M.id}{Cressie and Read (1985) M index}
  }
}
\source{
 Cressie, N (1991), \emph{Statistics for spatial data}. New York: Wiley, pp. 386--389; Cressie, N, Chan NH (1989) Spatial modelling of regional variables. \emph{Journal of the American Statistical Association}, 84, 393--401; Cressie, N, Read, TRC (1985) Do sudden infant deaths come in clusters? \emph{Statistics and Decisions} Supplement Issue 2, 333--349; \url{http://sal.agecon.uiuc.edu/datasets/sids.zip}.
}

\examples{
library(maptools)
sids.shp <- read.shape(system.file("shapes/sids.shp",
 package="maptools"))
sids <- sids.shp$att.data
sidspolys <- Map2poly(sids.shp, raw=FALSE)
sidscents <- get.Pcent(sids.shp)
data(nc.sids)
plot(sidspolys, border="grey", forcefill=FALSE)
plot(ncCR85.nb, sidscents, add=TRUE, col="blue")
plot(sidspolys, border="grey", forcefill=FALSE)
plot(ncCC89.nb, cbind(nc.sids$lon, nc.sids$lat),
 add=TRUE, col="blue")
}
\keyword{datasets}
