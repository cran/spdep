# Copyright 2001 by Roger Bivand
#


tri2nb <- function(coords, row.names = NULL) {
	require(tripack)
	n <- nrow(coords)
    	if (!is.null(row.names)) if(length(row.names) != n)
        	stop("row.names wrong length")
    	if (is.null(row.names)) row.names <- as.character(1:n)
	tri <- tri.mesh(x=coords[,1], y=coords[,2])
	nb <- neighbours(tri)
 	attr(nb, "region.id") <- row.names
	class(nb) <- "nb"
	attr(nb, "tri") <- TRUE
	attr(nb, "call") <- match.call()
	nb <- sym.attr.nb(nb)
	invisible(nb)
}

