# Copyright 2001 by Roger Bivand 
#


diffnb <- function(x, y, verbose=TRUE) {
	if (class(x) != "nb") stop("not a neighbours list")
	if (class(y) != "nb") stop("not a neighbours list")
	n <- length(x)
	if(n != length(y)) stop("lengths differ")
	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		xi <- x[[i]]
		yi <- y[[i]]
		xt <- xi %in% yi
		yt <- yi %in% xi
		if (!(all(xt) && all(yt))) {
			res[[i]] <- sort(unique(c(xi[which(!xt)],
				yi[which(!yt)])))
			if(verbose)
				cat("Neighbour difference for region:",
				i, "in relation to:", res[[i]], "\n")
		}
	}
	class(res) <- "nb"
	attr(res, "region.id") <- attr(x, "region.id")
	attr(res, "call") <- match.call()
	invisible(res)
}	
	
