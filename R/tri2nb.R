# Copyright 2001-4 by Roger Bivand
#


tri2nb <- function(coords, row.names = NULL) {
	require("tripack")
	n <- nrow(coords)
	if (n < 3) stop("too few coordinates")
#	left <- function(x) {
#		res <- (x[3,1]-x[2,1])*(x[1,2]-x[2,2]) >= 
#			(x[1,1]-x[2,1])*(x[3,2]-x[2,2])
#		res
#	}
#	if (left(coords[1:3,])) stop("first three coordinates collinear")
    	if (!is.null(row.names)) {
		if(length(row.names) != n)
            		stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names))
	    		stop("non-unique row.names given")
    	}
    	if (is.null(row.names)) row.names <- as.character(1:n)
	tri <- tri.mesh(x=coords[,1], y=coords[,2])
	nb <- neighbours(tri)
 	attr(nb, "region.id") <- row.names
	class(nb) <- "nb"
	attr(nb, "tri") <- TRUE
	attr(nb, "call") <- match.call()
	nb <- sym.attr.nb(nb)
	nb
}

