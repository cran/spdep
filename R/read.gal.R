# Copyright 2001 by Roger Bivand
#


read.gal <- function(file, row.names=NULL) 
{
	con <- file(file, open="r")
	n <- readLines(con, 1)
	n <- as.integer(n)
	if (n < 1) stop("Non-positive number of regions")
    	if (!is.null(row.names)) if(length(row.names) != n)
        	stop("row.names wrong length")
    	if (is.null(row.names)) row.names <- as.character(1:n)
	res <- vector(mode="list", length=n)
	for (i in 1:n) {
		line <- strsplit(readLines(con, 1), " ")
		x <- na.omit(as.integer(unlist(line)))
		if(x[1] != i) stop("GAL file corrupted")
		line <- strsplit(readLines(con, 1), " ")
		y <- na.omit(as.integer(unlist(line)))
		if(length(y) != x[2]) stop("GAL file corrupted")
		if(any(y < 0) || any(y > n))
			stop("GAL file corrupted")
		res[[i]] <- sort(y)
	}
	close(con)
	class(res) <- "nb"
    	attr(res, "region.id") <- row.names
	attr(res, "gal") <- TRUE
	attr(res, "call") <- TRUE
	res <- sym.attr.nb(res)
	invisible(res)
}

