# Copyright 2001-3 by Roger Bivand
#


nbdists <- function(nb, coords, lonlat=FALSE) {
	if (!inherits(nb, "nb")) 
        	stop("Not a neighbours list")
	if (!is.numeric(coords)) stop("Data non-numeric")
	if (!is.matrix(coords)) 
            stop("Data not in matrix form")
        if (any(is.na(coords))) 
            stop("Data include NAs")
	if (!is.double(coords)) storage.mode(coords) <- "double"
	n.nb <- length(nb)
	np <- nrow(coords)
        if (np != n.nb) 
            stop("Number of coords not equal to number of regions")
        dimension <- ncol(coords)
        dlist <- .Call("nbdists", nb, as.matrix(coords), as.integer(np), 
            as.integer(dimension), as.integer(lonlat), PACKAGE="spdep")
	attr(dlist[[1]], "call") <- match.call()
	invisible(dlist[[1]])
}

