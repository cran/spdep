# Copyright 2001 by Roger Bivand
#


summary.nb <- function(object, coords=NULL, ...) {
    nb <- object
    if (class(nb) != "nb") stop("Not a neighbours list")
    c.nb <- card(nb)
    n.nb <- length(nb)
    regids <- attr(nb, "region.id")
    if(is.null(regids)) regids <- as.character(1:n.nb)
    cat("Connectivity of", deparse(substitute(object)),
    	"with the following attributes:\n")
    print(str(attributes(nb)))
    cat("Number of regions:", n.nb, "\n")
    cat("Number of nonzero links:", sum(c.nb), "\n")
    cat("Percentage nonzero weights:", (100*sum(c.nb))/(n.nb*n.nb), "\n")
    cat("Average number of links:", mean(c.nb), "\n")
    cat("Link number distribution:\n")
    print(table(c.nb, deparse.level=0))
    if(any(c.nb == 0)) cat(length(c.nb[c.nb == 0]), " region", 
        ifelse(length(c.nb[c.nb == 0]) < 2, "", "s"), " with no links:\n",
	paste(regids[which(c.nb == 0)], collapse=" "), "\n", sep="")
    if(any(c.nb > 0)) {
        min.nb <- min(c.nb[c.nb > 0])
        cat(length(c.nb[c.nb == min.nb]), " least connected region",
	    ifelse(length(c.nb[c.nb == min.nb]) < 2, "", "s"), ":\n",
	    paste(regids[which(c.nb == min.nb)], collapse=" "), " with ",
	    min.nb, " link", ifelse(min.nb < 2, "", "s"), "\n", sep="")
        max.nb <- max(c.nb)
	cat(length(c.nb[c.nb == max.nb]), " most connected region",
	    ifelse(length(c.nb[c.nb == max.nb]) < 2, "", "s"), ":\n",
	    paste(regids[which(c.nb == max.nb)], collapse=" "), " with ",
	    max.nb, " link", ifelse(max.nb < 2, "", "s"), "\n", sep="")
    }
    if(!is.null(coords)) {
        if (!is.matrix(coords)) stop("Data not in matrix form")
        if (any(is.na(coords))) stop("Data include NAs")
        np <- nrow(coords)
	if(np != n.nb) stop("Number of coords not equal to number of regions")
        dimension <- ncol(coords)
	dlist <- .Call("nbdists", nb, as.matrix(coords), as.integer(np), 
	    as.integer(dimension), PACKAGE="spdep")[[1]]
	cat("Summary of link distances:\n")
	print(summary(unlist(dlist)))
	stem(unlist(dlist))
    }
}
