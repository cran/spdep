# Copyright 2004 by Roger Bivand 
#

nb2blocknb <- function(nb, ID, row.names = NULL) {
	if (!inherits(nb, "nb")) stop("not an nb object")
	nbNames <- as.character(attr(nb, "region.id"))
	entNames <- as.character(ID)
	if (!identical(sort(nbNames), sort(unique(entNames))))
		stop("names do not match exactly")
	n <- length(entNames)
	if (n < 1) stop("non-positive number of entities")
	if (!is.null(row.names)) {
		if (length(row.names) != n) 
			stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names)) 
		stop("non-unique row.names given")
	} else {
		row.names <- as.character(1:n)
	}
	inter <- lapply(as.list(nbNames), 
		function(x) which(match(entNames, x) == 1))

	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		ii <- match(entNames[i], nbNames)
		blocks <- c(ii, nb[[ii]])
		vec <- sort(unlist(inter[blocks]))
		res[[i]] <- vec[vec != i]
	}

	attr(res, "region.id") <- row.names
	class(res) <- "nb"
	attr(res, "block") <- TRUE
	attr(res, "call") <- match.call()
	res <- sym.attr.nb(res)
	invisible(res)
}


