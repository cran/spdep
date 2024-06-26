% Copyright 2001-2019 by Nicholas Lewin-Koh and Roger Bivand
\name{graphneigh}
\alias{gabrielneigh}
\alias{relativeneigh}
\alias{soi.graph}
\alias{plot.Gabriel}
\alias{plot.relative}
\alias{graph2nb}


\title{Graph based spatial weights}
\description{
Functions return a graph object containing a list with the vertex
coordinates and the to and from indices defining the edges. Some/all of these functions assume that the coordinates are not exactly regularly spaced. The helper
function \code{graph2nb} converts a graph
object into a neighbour list. The plot functions plot the graph objects.
}
\usage{
gabrielneigh(coords, nnmult=3)
relativeneigh(coords, nnmult=3)
%beta.skel(coords,beta)
soi.graph(tri.nb, coords, quadsegs=10)
graph2nb(gob, row.names=NULL,sym=FALSE)
\method{plot}{Gabriel}(x, show.points=FALSE, add=FALSE, linecol=par(col), ...)
\method{plot}{relative}(x, show.points=FALSE, add=FALSE, linecol=par(col),...)
}

\arguments{
  \item{coords}{matrix of region point coordinates or SpatialPoints object or \code{sfc} points object}
  \item{nnmult}{scaling factor for memory allocation, default 3; if higher values are required, the function will exit with an error; example below thanks to Dan Putler}
  \item{tri.nb}{a neighbor list created from tri2nb}
  \item{quadsegs}{number of line segments making a quarter circle buffer, see the \code{nQuadSegs} argument in \code{\link[sf]{geos_unary}}}
%  \item{beta}{the parameter for a beta skeleton}
  \item{gob}{a graph object created from any of the graph funtions}
  \item{row.names}{character vector of region ids to be added to the
    neighbours list as attribute \code{region.id}, default \code{seq(1,
      nrow(x))}}
  \item{sym}{a logical argument indicating whether or not neighbors
    should be symetric (if i->j then j->i)}
  \item{x}{object to be plotted}
  \item{show.points}{(logical) add points to plot}
  \item{add}{(logical) add to existing plot}
  \item{linecol}{edge plotting colour}
  \item{...}{further graphical parameters as in \code{par(..)}}
}
\details{
The graph functions produce graphs on a 2d point set that
%except for
%some values of \eqn{ \beta (\beta < 1)}{\beta (\beta < 1)} in the
%beta-skeleton, 
are all subgraphs of the Delaunay triangulation. The
relative neighbor graph is defined by the relation, x and y are neighbors if

\deqn{d(x,y) \le min(max(d(x,z),d(y,z))| z \in S)}{d(x,y) <= min(max(d(x,z),d(y,z))| z in S)}

where d() is the distance, S is the set of points and z is an arbitrary
point in S. The Gabriel graph is a subgraph of the delaunay
triangulation and has the relative neighbor graph as a sub-graph. The
relative neighbor graph is defined by the relation x and y are Gabriel
neighbors if

\deqn{d(x,y) \le min((d(x,z)^2 + d(y,z)^2)^{1/2} |z \in S)}{d(x,y) <= min((d(x,z)^2 + d(y,z)^2)^1/2 |z in S)}

where x,y,z and S are as before. The sphere of influence graph is
defined for a finite point set S, let \eqn{r_x} be the distance from point x
to its nearest neighbor in S, and \eqn{C_x} is the circle centered on x. Then
x and y are SOI neigbors iff \eqn{C_x} and \eqn{C_y} intersect in at
least 2 places. From 2016-05-31, Computational Geometry in C code replaced by calls to functions in \pkg{dbscan} and \pkg{sf}; with a large \code{quadsegs=} argument, the behaviour of the function is the same, otherwise buffer intersections only closely approximate the original function.

%The \eqn{beta}{\beta} 
See \code{\link{card}} for details of \dQuote{nb} objects.
}
\value{
A list of class \code{Graph} with the following elements
  \item{np}{number of input points}
  \item{from}{array of origin ids}
  \item{to}{array of destination ids}
  \item{nedges}{number of edges in graph}
  \item{x}{input x coordinates}
  \item{y}{input y coordinates}
The helper functions return an \code{nb} object with a list of integer
 vectors containing neighbour region number ids.
}

\references{
  Matula, D. W. and Sokal R. R. 1980, Properties of Gabriel
  graphs relevant to geographic variation research and the clustering of
  points in the plane, Geographic Analysis, 12(3), pp. 205-222.

  Toussaint, G. T. 1980, The relative neighborhood graph of a finite
  planar set, Pattern Recognition, 12(4), pp. 261-268.

  Kirkpatrick, D. G. and Radke, J. D. 1985, A framework for
  computational morphology. In Computational Geometry,
  Ed. G. T. Toussaint, North Holland.

  
  }

\author{Nicholas Lewin-Koh \email{nikko@hailmail.net}}

\seealso{\code{\link{knearneigh}}, \code{\link{dnearneigh}},
\code{\link{knn2nb}}, \code{\link{card}}}

\examples{
columbus <- st_read(system.file("shapes/columbus.gpkg", package="spData")[1], quiet=TRUE)
sf_obj <- st_centroid(st_geometry(columbus), of_largest_polygon)
sp_obj <- as(sf_obj, "Spatial")
coords <- st_coordinates(sf_obj)
suppressMessages(col.tri.nb <- tri2nb(coords))
col.gab.nb <- graph2nb(gabrielneigh(coords), sym=TRUE)
col.rel.nb <- graph2nb(relativeneigh(coords), sym=TRUE)
par(mfrow=c(2,2))
plot(st_geometry(columbus), border="grey")
plot(col.tri.nb,coords,add=TRUE)
title(main="Delaunay Triangulation", cex.main=0.6)
plot(st_geometry(columbus), border="grey")
plot(col.gab.nb, coords, add=TRUE)
title(main="Gabriel Graph", cex.main=0.6)
plot(st_geometry(columbus), border="grey")
plot(col.rel.nb, coords, add=TRUE)
title(main="Relative Neighbor Graph", cex.main=0.6)
plot(st_geometry(columbus), border="grey")
if (require("dbscan", quietly=TRUE)) {
  col.soi.nb <- graph2nb(soi.graph(col.tri.nb,coords), sym=TRUE)
  plot(col.soi.nb, coords, add=TRUE)
  title(main="Sphere of Influence Graph", cex.main=0.6)
}
par(mfrow=c(1,1))
col.tri.nb_sf <- tri2nb(sf_obj)
all.equal(col.tri.nb, col.tri.nb_sf, check.attributes=FALSE)
col.tri.nb_sp <- tri2nb(sp_obj)
all.equal(col.tri.nb, col.tri.nb_sp, check.attributes=FALSE)
if (require("dbscan", quietly=TRUE)) {
  col.soi.nb_sf <- graph2nb(soi.graph(col.tri.nb, sf_obj), sym=TRUE)
  all.equal(col.soi.nb, col.soi.nb_sf, check.attributes=FALSE)
  col.soi.nb_sp <- graph2nb(soi.graph(col.tri.nb, sp_obj), sym=TRUE)
  all.equal(col.soi.nb, col.soi.nb_sp, check.attributes=FALSE)
}
col.gab.nb_sp <- graph2nb(gabrielneigh(sp_obj), sym=TRUE)
all.equal(col.gab.nb, col.gab.nb_sp, check.attributes=FALSE)
col.gab.nb_sf <- graph2nb(gabrielneigh(sf_obj), sym=TRUE)
all.equal(col.gab.nb, col.gab.nb_sf, check.attributes=FALSE)
col.rel.nb_sp <- graph2nb(relativeneigh(sp_obj), sym=TRUE)
all.equal(col.rel.nb, col.rel.nb_sp, check.attributes=FALSE)
col.rel.nb_sf <- graph2nb(relativeneigh(sf_obj), sym=TRUE)
all.equal(col.rel.nb, col.rel.nb_sf, check.attributes=FALSE)
dx <- rep(0.25*0:4,5)
dy <- c(rep(0,5),rep(0.25,5),rep(0.5,5), rep(0.75,5),rep(1,5))
m <- cbind(c(dx, dx, 3+dx, 3+dx), c(dy, 3+dy, dy, 3+dy))
cat(try(res <- gabrielneigh(m), silent=TRUE), "\n")
res <- gabrielneigh(m, nnmult=4)
summary(graph2nb(res))
grd <- as.matrix(expand.grid(x=1:5, y=1:5)) #gridded data
r2 <- gabrielneigh(grd)
set.seed(1)
grd1 <- as.matrix(expand.grid(x=1:5, y=1:5)) + matrix(runif(50, .0001, .0006), nrow=25)
r3 <- gabrielneigh(grd1)
opar <- par(mfrow=c(1,2))
plot(r2, show=TRUE, linecol=2)
plot(r3, show=TRUE, linecol=2)
par(opar)
# example of reading points with readr::read_csv() yielding a tibble
load(system.file("etc/misc/coords.rda", package="spdep"))
class(coords)
graph2nb(gabrielneigh(coords))
graph2nb(relativeneigh(coords))
}
\keyword{spatial}
