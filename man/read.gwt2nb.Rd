% Copyright 2003 by Roger S. Bivand
\name{read.gwt2nb}
\alias{read.gwt2nb}
\alias{write.sn2gwt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Read and write GeoDa GWT spatial neighbour files}
\description{
  
}
\usage{
read.gwt2nb(file, region.id=NULL)
write.sn2gwt(sn, file, shpfile=NULL, ind=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{name of file with GWT lattice data}
  \item{region.id}{region IDs}
  \item{sn}{a \code{spatial.neighbour} object}
  \item{shpfile}{Shapefile name taken from GWT file for this dataset}
  \item{ind}{region id indicator variable name}
}
\details{
  
}
\value{
  
}

\references{Luc Anselin (2003) \emph{GeoDa 0.9 User's Guide}, pp. 80--81, Spatial Analysis Laboratory, Department of Agricultural and Consumer Economics, University of Illinois, Urbana-Champaign, \url{http://agec221.agecon.uiuc.edu/csiss/pdf/geoda093.pdf}}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
 
\seealso{\code{\link{read.gal}}}

\examples{
data(baltimore)
STATION <- baltimore$STATION
gwt1 <- read.gwt2nb(system.file("etc/weights/baltk4.GWT", package="spdep")[1],
 STATION)
listw1 <- nb2listw(gwt1, style="B", glist=attr(gwt1, "GeoDa")$dist)
tmpGWT <- tempfile()
write.sn2gwt(listw2sn(listw1), tmpGWT)
gwt2 <- read.gwt2nb(tmpGWT, STATION)
diffnb(gwt1, gwt2)
}
\keyword{spatial}
