# Copyright 2002-3 by Roger Bivand 
#

eigenw <- function(listw, quiet=TRUE)
{
	if(!inherits(listw, "listw")) stop("not a listw object")
	w <- listw2mat(listw)
	e <- eigen(w, only.values=TRUE)$values
	if (!quiet) {
		cat("Largest eigenvalue:", 
		if(is.complex(e)) max(Re(e)) else max(e),
		"Sum of eigenvalues:", sum(e), "\n")
	}
	e
}

