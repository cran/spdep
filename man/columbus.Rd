% Copyright 2001 by Roger S. Bivand
\name{columbus}
\alias{columbus}
\alias{bbs}
\alias{col.gal.nb}
\alias{coords}
\alias{polys}
\non_function{}
\title{Columbus OH spatial analysis data set}
\usage{data(columbus)}
\description{
The \code{columbus} data frame has 49 rows and 22 columns.
Unit of analysis: 49 neighbourhoods in Columbus, OH, 1980 data.
In addition the data set includes a \code{polylist} object \code{polys}
with the boundaries of the neighbourhoods, a matrix of polygon centroids
\code{coords}, a matrix \code{bbs} with
the bounding boxes of the polygons, and \code{col.gal.nb}, the neighbours list
from read.gal("spweights/data/col-queen.gal").
}
\format{
  This data frame contains the following columns:
  \describe{
    \item{AREA}{computed by ArcView}
    \item{PERIMETER}{computed by ArcView}
    \item{COLUMBUS.}{internal polygon ID (ignore)}
    \item{COLUMBUS.I}{another internal polygon ID (ignore)}
    \item{POLYID}{yet another polygon ID}
    \item{NEIG}{neighborhood id value (1-49);
             conforms to id value used in Spatial Econometrics book.}
    \item{HOVAL}{housing value (in \$1,000)}
    \item{INC}{household income (in \$1,000)}
    \item{CRIME}{residential burglaries and vehicle thefts per thousand
        households in the neighborhood}
    \item{OPEN}{open space in neighborhood}
    \item{PLUMB}{percentage housing units without plumbin}
    \item{DISCBD}{distance to CBD}
    \item{X}{x coordinate (in arbitrary digitizing units, not polygon coordinates)}
    \item{Y}{y coordinate (in arbitrary digitizing units, not polygon coordinates)}
    \item{AREA}{neighborhood area (computed by SpaceStat)}
    \item{NSA}{north-south dummy (North=1)}
    \item{NSB}{north-south dummy (North=1)}
    \item{EW}{east-west dummy (East=1)}
    \item{CP}{core-periphery dummy (Core=1)}
    \item{THOUS}{constant=1,000}
    \item{NEIGNO}{NEIG+1,000, alternative neighborhood id value}
    \item{PERIM}{polygon perimeter (computed by SpaceStat)}
  }
}
\details{
The row names of \code{columbus} and the \code{region.id} attribute of
\code{polys} are set to \code{columbus$NEIGNO}.
}
\source{
Anselin, Luc.  1988.  Spatial econometrics: methods and
models.  Dordrecht: Kluwer Academic, Table 12.1 p. 189.
\url{http://www.spacestat.com},
\url{ftp://cho.pol.uiuc.edu/pub/spacestat/data/columbus.ZIP}
}
\note{
All source data files prepared by
Luc Anselin, Bruton Center for Development Studies
University of Texas at Dallas
P.O. Box 830688 - GR 31
Richardson, TX 75083-0688
e-mail: lanselin@utdallas.edu
}

\examples{
data(columbus)
crime.f <- as.ordered(cut(columbus$CRIME,
  breaks=quantile(columbus$CRIME, probs=seq(0,1,0.2)),
  include.lowest=TRUE))
colours <- c("salmon1", "salmon2", "red3", "brown", "black")
plot(bbs[,1], bbs[,4], xlab="", ylab="", asp=1, type="n",
  xlim=range(c(bbs[,1],bbs[,3])), ylim=range(c(bbs[,2],bbs[,4])))
for (i in 1:length(polys)) 
  polygon(polys[[i]], col=colours[codes(crime.f[i])])
legend(x=c(5.8, 7.1), y=c(13, 14.5), legend=levels(crime.f),
  fill=colours, bty="n")
title(main=paste("Columbus OH: residential burglaries and vehicle\nthefts",
  "per thousand households, 1980"))
}
\keyword{datasets}
