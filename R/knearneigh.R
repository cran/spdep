# Copyright 2001-2010 by Roger S. Bivand. 
# Upgrade to sp classes February 2007
# Added RANN April 2010

knearneigh <- function(x, k=1, longlat=NULL, RANN=TRUE)
{
    if (inherits(x, "SpatialPoints")) {
        if ((is.null(longlat) || !is.logical(longlat)) 
 	   && !is.na(is.projected(x)) && !is.projected(x)) {
           longlat <- TRUE
        } else longlat <- FALSE
        x <- coordinates(x)[, 1:2]
    } else if (is.null(longlat) || !is.logical(longlat)) longlat <- FALSE
    if (!is.numeric(x)) stop("Data non-numeric")
    if (!is.matrix(x)) stop("Data not in matrix form")
    stopifnot(ncol(x) == 2)
    if (any(is.na(x))) stop("Data include NAs")
    if (longlat) {
        bb <- bbox(x)
        if (!sp:::.ll_sanity(bb))
            warning("Coordinates are not geographical: longlat argument wrong")
    }
    if (!is.double(x)) storage.mode(x) <- "double"
    np <- nrow(x)
    dimension <- ncol(x)
    if (dimension != 2) stop("Only 2D data accepted")
    if (k >= np) stop("Fewer data points than k")
    if (RANN && !longlat && require(RANN, quietly=TRUE)) {
        xx <- cbind(x, out=rep(0, nrow(x)))
        out <- as.matrix(nn(xx, p=k)$nn.idx)
        dimnames(out) <- NULL
        res <- list(nn=out, np=np, k=k, dimension=dimension, x=x)
    } else {
        xx <- c(x[,1], x[,2])
        storage.mode(xx) <- "double"
        nn <- integer(np*k)
        dnn <- double(np*k)
        z <- .C("knearneigh", k=as.integer(k), np=as.integer(np),
            dimension=as.integer(dimension),
            xx=xx, nn=as.integer(nn), dnn=dnn,
	    as.integer(longlat), PACKAGE="spdep")
        res <- list(nn=matrix(z$nn, np, k, byrow=TRUE), np=np, k=k,
    	    dimension=dimension, x=x)
    }
    class(res) <- "knn"
    attr(res, "call") <- match.call()
    res
}
