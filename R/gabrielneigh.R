# Copyright 2001 by Nicholas Lewin-Koh 
#

gabrielneigh <- function(coords) {
    x <- coords
    if (!is.matrix(x)) stop("Data not in matrix form")
    if (any(is.na(x))) stop("Data cannot include NAs")
    np <- nrow(x)
    if(ncol(x)!=2) stop("Planar graphs only work in 2d")
    g1<-g2<-rep(0,np*3)
    nogab <- 0
    z <- .C("compute_gabriel", np=as.integer(np), from=as.integer(g1),
             to=as.integer(g2), nedges=as.integer(nogab),
             x=as.double(x[,1]), y=as.double(x[,2]), PACKAGE="spdep")
    z$from<-z$from[1:z$nedges]
    z$to<-z$to[1:z$nedges]
    attr(z, "call") <- match.call()
    class(z)<-c("Graph","Gabriel")
    invisible(z)
}

plot.Gabriel<-function(x, show.points=FALSE, add=FALSE,
                       linecol=par(col), ...)
{
  if(!add) plot(gab$x,gab$y,type='n',...)
  segments(gab$x[gab$from], gab$y[gab$from],
           gab$x[gab$to], gab$y[gab$to], col=linecol)
  if(show.points) points(gab$x,gab$y)
}

