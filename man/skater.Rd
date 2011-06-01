\name{skater}
\alias{skater}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial 'K'luster Analysis by Tree Edge Removal}
\description{
  This function implements a SKATER procedure for spatial clustering
  analysis. This procedure essentialy begins with an edges set, a data set and
  a number of cuts. The output is an object of 'skater' class and is
  valid for input again.
}
\usage{
skater(edges, data, ncuts, crit, vec.crit, method = c("euclidean", 
    "maximum", "manhattan", "canberra", "binary", "minkowski", 
    "mahalanobis"), p = 2, cov, inverted = FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{edges}{A matrix with 2 colums with each row is an edge}
  \item{data}{A data.frame with data observed over nodes.}
  \item{ncuts}{The number of cuts}
  \item{crit}{A numeric or integer with criteria for groups.
    Example: minimum population size.}
  \item{vec.crit}{A vector for evaluating critera.}
  \item{method}{Character or function to declare distance method.
    If \code{method} is character, method must be "mahalanobis" or
    "euclidean", "maximum", "manhattan", "canberra", "binary"
    or "minkowisk".
    If \code{method} is one of "euclidean", "maximum",
    "manhattan", "canberra", "binary" or "minkowisk", see
    \code{\link{dist}} for details,
     because this function as used to compute the distance.
     If \code{method="mahalanobis"}, the mahalanobis distance
     is computed between neighbour areas.
     If \code{method} is a \code{function}, this function is
     used to compute the distance.}
  \item{p}{The power of the Minkowski distance.}
  \item{cov}{The covariance matrix used to compute the mahalanobis
   distance.}
  \item{inverted}{logical.  If 'TRUE', 'cov' is supposed to contain the
          inverse of the covariance matrix.}
}
\details{
  to do
}
\value{
  A object of \code{skater} class with:
  
  \item{groups}{A vector with length equal the number of nodes. Each
    position identifies the group of node} 
  \item{edges.groups}{A list of length equal the number of groups with
    each element is a set of edges} 
  \item{not.prune}{A vector identifying the groups with are not
    candidates to partition.}
  \item{candidates}{A vector identifying the groups with are candidates
    to partition.} 
  \item{ssto}{The total dissimilarity in each step of edge removal.}

}
\references{
  Assuncao, R.M., Lage J.P., and Reis, E.A. (2002). Analise de
  conglomerados espaciais via arvore geradora minima. Revista Brasileira
  de Estatistica, 62, 1-23. 

  Assuncao, R. M, Neves, M. C., Camara, G. and Freitas, C. da
  C. (2006). Efficient regionalization techniques for socio-economic
  geographical units using minimum spanning trees. International Journal
  of Geographical Information Science Vol. 20, No. 7, August 2006,
  797-811 
}
\author{Renato M. Assuncao and Elias T. Krainski}
%%\note{}
\seealso{See Also as \code{\link{mstree}}}
\examples{
### loading data
bh <- readShapePoly(system.file("etc/shapes/bhicv.shp",
      package="spdep")[1])
### data standardized 
dpad <- data.frame(scale(bh@data[,5:8]))

### neighboorhod list
bh.nb <- poly2nb(bh)

### calculating costs
lcosts <- nbcosts(bh.nb, dpad)

### making listw
nb.w <- nb2listw(bh.nb, lcosts, style="B")

### find a minimum spanning tree
mst.bh <- mstree(nb.w,5)

### the mstree plot
par(mar=c(0,0,0,0))
plot(mst.bh, coordinates(bh), col=2,       
     cex.lab=.7, cex.circles=0.035, fg="blue")
plot(bh, border=gray(.5), add=TRUE)

### three groups with no restriction
res1 <- skater(mst.bh[,1:2], dpad, 2)

### thee groups with minimum population 
res2 <- skater(mst.bh[,1:2], dpad, 2, 200000, bh@data$Pop)

### thee groups with minimun number of areas
res3 <- skater(mst.bh[,1:2], dpad, 2, 3, rep(1,nrow(bh@data)))

### groups frequency
table(res1$groups)
table(res2$groups)
table(res3$groups)

### the skater plot
par(mar=c(0,0,0,0))
plot(res1, coordinates(bh), cex.circles=0.035, cex.lab=.7)

### more one partition
res1b <- skater(res1, dpad, 1)

### length groups frequency
table(res1$groups)
table(res1b$groups)

### the skater plot, using other colors
plot(res1b, coordinates(bh), cex.circles=0.035, cex.lab=.7,
     groups.colors=colors()[(1:length(res1b$ed))*10])

### the Spatial Polygons plot
plot(bh, col=heat.colors(4)[res1b$groups])

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
\keyword{tree}% __ONLY ONE__ keyword per line
