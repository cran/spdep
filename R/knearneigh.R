# Copyright 2001-4 by Roger S. Bivand. 
#

knearneigh <- function(x, k=1, lonlat=FALSE)
{
    if (!is.numeric(x)) stop("Data non-numeric")
    if (!is.matrix(x)) stop("Data not in matrix form")
    if (any(is.na(x))) stop("Data include NAs")
    if (!is.double(x)) storage.mode(x) <- "double"
    np <- nrow(x)
    dimension <- ncol(x)
    if (dimension != 2) stop("Only 2D data accepted")
    if (k >= np) stop("Fewer data points than k")
    xx <- c(x[,1], x[,2])
    nn <- integer(np*k)
    dnn <- double(np*k)
    z <- .C("knearneigh", k=as.integer(k), np=as.integer(np),
        dimension=as.integer(dimension),
        xx=as.double(xx), nn=as.integer(nn), dnn=as.double(dnn),
	as.integer(lonlat), PACKAGE="spdep")
    res <- list(nn=matrix(z$nn, np, k, byrow=TRUE), np=np, k=k,
    	dimension=dimension, x=x)
    class(res) <- "knn"
    attr(res, "call") <- match.call()
    invisible(res)
}
