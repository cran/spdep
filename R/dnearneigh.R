# Copyright 2000-3 by Roger S. Bivand. 
#

dnearneigh <- function(x, d1, d2, row.names=NULL, lonlat=FALSE) {
    if (!is.numeric(x)) stop("Data non-numeric")
    if (!is.matrix(x)) stop("Data not in matrix form")
    if (any(is.na(x))) stop("Data include NAs")
    if (!is.double(x)) storage.mode(x) <- "double"
    np <- nrow(x)
    if (!is.null(row.names)) {
	if(length(row.names) != np)
            stop("row.names wrong length")
	if (length(unique(row.names)) != length(row.names))
	    stop("non-unique row.names given")
    }
    if (is.null(row.names)) row.names <- as.character(1:np)
    dimension <- ncol(x)
    if (dimension > 2) stop("Only 2D data accepted")
    md <- 0
    if (d1 < 0) d1 <- 0.0
    if (!lonlat) {
	for (i in 1:dimension) md <- sum(md, (diff(range(x[,i]))^2))
    	if (d2 > sqrt(md)) d2 <- sqrt(md)
    }
    z <- .Call("dnearneigh", as.double(d1), as.double(d2), as.integer(np),
        as.integer(dimension), as.double(x), as.integer(lonlat), 
	PACKAGE="spdep")
    attr(z[[1]], "region.id") <- row.names
    attr(z[[1]], "call") <- match.call()
    attr(z[[1]], "dnn") <- c(d1, d2)
    z[[1]] <- sym.attr.nb(z[[1]])
    invisible(z[[1]])
}
