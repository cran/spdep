% Copyright 2003 by Roger S. Bivand
\name{Map2poly}
\alias{Map2poly}
\alias{Map2bbs}
\alias{shape2poly}
\alias{shape2bbs}
\alias{convert.pl}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Create polygon lists and bounding boxes from imported shapefiles}
\description{
  Map2poly() is a simple function to make imported GIS vector polygons into a "polylist" object from a "Map" object, and Map2bbs() retrieves polygon bounding boxes; analogous functions for the "shapefiles" package are shape2poly() and shape2bbs() (thanks to Stéphane Dray for his contribution); convert.pl() serves to convert the deprecated "multipart" "polylist" format to an NA-separated format now used by plotpolys() and poly2nb().
}
\usage{
Map2poly(Map, region.id = NULL)
Map2bbs(Map)
shape2poly(shape, region.id = NULL)
shape2bbs(shape)
convert.pl(pl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Map}{object of class \code{"Map"}, imported from a GIS polygon vector file using unpublished package "maptools"}
  \item{shape}{list returned by \code{read.shapefile()}, imported from a GIS polygon vector file using the package "shapefiles"}
  \item{region.id}{character vector of region ids to be added to the neighbours list as attribute \code{region.id}}
  \item{pl}{list of old-style multipart polygons of class \code{polylist}}
}
\details{
  \url{http://sal.agecon.uiuc.edu/datasets/columbus.zip}
}

\value{
  \code{Map2poly}, \code{shape2poly}, \code{convert.pl} return polylist objects, lists of polygon boundary coordinates (divided by NA where the polygon object is represented by more than one polygon;
  \code{Map2bbs}, \code{shape2bbs} return bounding box matrixes, c(x1, y1, x2, y2) with coordinates for the lower left corner and upper right corner.
}

\author{Roger Bivand and Stéphane Dray, \email{Roger.Bivand@nhh.no}}
 
\seealso{\code{\link{read.shape}}, \code{\link{read.shapefile}}}


\examples{
shp <- try(library(shapefiles))
if (class(shp) != "try-error") {
  try1 <- read.shapefile(paste(system.file("etc/columbus",
    package="spdep")[1], "columbus", sep="/"))
  shppolys <- shape2poly(try1, as.character(try1$dbf$dbf$NEIGNO))
  shpbbs <- shape2bbs(try1)
  plotpolys(shppolys, shpbbs)
  title(main="Polygons for Columbus OH from shapefiles package")
}
map <- try(library(maptools))
if (class(map) != "try-error") {
  try2 <- read.shape(system.file("etc/columbus/columbus.shp", 
    package="spdep")[1])
  mappolys <- Map2poly(try2, as.character(try2$att.data$NEIGNO))
  mapbbs <- Map2bbs(try2)
  plotpolys(mappolys, mapbbs)
  title(main="Polygons for Columbus OH from maptools package")
  if (class(shp) != "try-error") {
    plotpolys(shppolys, shpbbs)
    plotpolys(mappolys, mapbbs, add=TRUE, border="red", lty="dotted")
    title(main="Polygons for Columbus OH overplotted from both packages")
  }
}
}

}
\keyword{spatial}
