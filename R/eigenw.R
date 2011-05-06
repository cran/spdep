# Copyright 2002-11 by Roger Bivand 
#

eigenw <- function(listw, quiet=NULL)
{
	if(!inherits(listw, "listw")) stop("not a listw object")
        if (is.null(quiet)) quiet <- !get("verbose", env = .spdepOptions)
        stopifnot(is.logical(quiet))

	w <- listw2mat(listw)
	e <- eigen(w, only.values=TRUE)$values
	if (!quiet) {
		cat("Largest eigenvalue:", 
# modified 110414 RSB
		if(is.complex(e)) {max(Re(e[which(Im(e) == 0)]))} else max(e),
		"Sum of eigenvalues:", sum(e), "\n")
	}
	e
}

