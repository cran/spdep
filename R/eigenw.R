# Copyright 2002 by Roger Bivand 
#

eigenw <- function(listw, quiet=TRUE)
{
	if(class(listw) != "listw") stop("not a listw object")
	w <- listw2mat(listw)
	e <- eigen(w, only.values=TRUE)$values
	if (is.complex(e)) e <- Re(e)
	if (!quiet) {
		cat("Largest eigenvalue:", max(e),
		"Sum of eigenvalues:", sum(e), "\n")
	}
	e
}

