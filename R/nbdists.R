# Copyright 2001 by Roger Bivand
#


nbdists <- function(nb, coords) {
	if (class(nb) != "nb") 
        	stop("Not a neighbours list")
	if (!is.matrix(coords)) 
            stop("Data not in matrix form")
        if (any(is.na(coords))) 
            stop("Data include NAs")
	n.nb <- length(nb)
	np <- nrow(coords)
        if (np != n.nb) 
            stop("Number of coords not equal to number of regions")
        dimension <- ncol(coords)
        dlist <- .Call("nbdists", nb, as.matrix(coords), as.integer(np), 
            as.integer(dimension), PACKAGE="spdep")
	attr(dlist[[1]], "call") <- match.call()
	invisible(dlist[[1]])
}

