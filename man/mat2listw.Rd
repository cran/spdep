% Copyright 2003-2010 by Roger S. Bivand
\name{mat2listw}
\alias{mat2listw}

\title{Convert a square spatial weights matrix to a weights list object}
\description{
  The function converts a square spatial weights matrix, 
optionally a sparse matrix to a weights list 
object, optionally adding region IDs from the row names of the matrix, as a 
sequence of numbers 1:nrow(x), or as given as an argument. The style can be imposed by rebuilting the weights list object internally.
}
\usage{
mat2listw(x, row.names = NULL, style=NULL, zero.policy = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A square non-negative matrix with no NAs representing spatial 
weights; may be a matrix of class \dQuote{sparseMatrix}}
  \item{row.names}{row names to use for region IDs}
  \item{style}{default NULL, missing, set to "M" and warning given; if not "M", passed to \code{\link{nb2listw}} to re-build the object}
  \item{zero.policy}{default NULL, use global option value; if FALSE stop with error for any empty neighbour sets, if TRUE permit the weights list to be formed with zero-length weights vectors}
}

\value{
A \code{listw} object with the following members:
  \item{style}{"M", meaning matrix style, underlying style unknown, or assigned style argument in rebuilt object}
  \item{neighbours}{the derived neighbours list}
  \item{weights}{the weights for the neighbours derived from the matrix}
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

 
\seealso{\code{\link{nb2listw}}, \code{\link{nb2mat}}}
\examples{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
col005 <- dnearneigh(st_coordinates(st_centroid(st_geometry(columbus),
 of_largest_polygon=TRUE)), 0, 0.5, as.character(columbus$NEIGNO))
summary(col005)
col005.w.mat <- nb2mat(col005, style="W", zero.policy=TRUE)
try(col005.w.b <- mat2listw(col005.w.mat, style="W"))
col005.w.b <- mat2listw(col005.w.mat, style="W", zero.policy=TRUE)
summary(col005.w.b$neighbours)
diffnb(col005, col005.w.b$neighbours)
col005.w.mat.3T <- kronecker(diag(3), col005.w.mat)
col005.w.b.3T <- mat2listw(col005.w.mat.3T, style="W", zero.policy=TRUE)
summary(col005.w.b.3T$neighbours)
run <- FALSE
if (require("spatialreg", quiet=TRUE)) run <- TRUE
if (run) {
W <- as(nb2listw(col005, style="W", zero.policy=TRUE), "CsparseMatrix")
try(col005.spM <- mat2listw(W))
col005.spM <- mat2listw(W, style="W", zero.policy=TRUE)
summary(col005.spM$neighbours)
}
if (run) {
diffnb(col005, col005.spM$neighbours)
}
if (run && require("Matrix", quiet=TRUE)) {
IW <- kronecker(Diagonal(3), W)
col005.spM.3T <- mat2listw(as(IW, "CsparseMatrix"), style="W", zero.policy=TRUE)
summary(col005.spM.3T$neighbours)
}
}
\keyword{spatial}
