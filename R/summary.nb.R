# Copyright 2001-3 by Roger Bivand
#


summary.nb <- function(object, coords=NULL, lonlat=FALSE, scale=1, ...) {
    nb <- object
    if (!inherits(nb, "nb")) stop("Not a neighbours list")
    c.nb <- card(nb)
    n.nb <- length(nb)
    regids <- attr(nb, "region.id")
    if(is.null(regids)) regids <- as.character(1:n.nb)
    print.nb(object)
    cat("Link number distribution:\n")
    print(table(c.nb, deparse.level=0))
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
	    as.integer(dimension), as.integer(lonlat), PACKAGE="spdep")[[1]]
	cat("Summary of link distances:\n")
	print(summary(unlist(dlist)))
	stem(unlist(dlist), scale=scale)
    }
}


print.nb <- function(x, ...) {
    nb <- x
    if (!inherits(nb, "nb")) stop("Not a neighbours list")
    c.nb <- card(nb)
    n.nb <- length(nb)
    regids <- attr(nb, "region.id")
    if(is.null(regids)) regids <- as.character(1:n.nb)
    cat("Neighbour list object:\n")
    cat("Number of regions:", n.nb, "\n")
    cat("Number of nonzero links:", sum(c.nb), "\n")
    cat("Percentage nonzero weights:", (100*sum(c.nb))/(n.nb^2), "\n")
    cat("Average number of links:", mean(c.nb), "\n")
    if(any(c.nb == 0)) cat(length(c.nb[c.nb == 0]), " region", 
        ifelse(length(c.nb[c.nb == 0]) < 2, "", "s"), " with no links:\n",
	paste(regids[which(c.nb == 0)], collapse=" "), "\n", sep="")
    res <- is.symmetric.nb(x, verbose=FALSE)
    if (!res) cat("Non-symmetric neighbours list\n")
    invisible(x)
}

summary.listw <- function(object, coords=NULL, lonlat=FALSE, 
	zero.policy=FALSE, scale=1, ...) {
	cat("Characteristics of weights list object:\n")
	summary(object$neighbours, coords=coords, lonlat=lonlat, 
		scale=scale, ...)
	style <- object$style
	cat(paste("\nWeights style:", style, "\n"))
	cat("Weights constants summary:\n")
	print(data.frame(rbind(unlist(spweights.constants(object,
		zero.policy=zero.policy))[c(1, 5:8)]), row.names=style))

}

print.listw <- function(x, zero.policy=FALSE, ...) {
	cat("Characteristics of weights list object:\n")
	print.nb(x$neighbours, ...)
	style <- x$style
	cat(paste("\nWeights style:", style, "\n"))
	cat("Weights constants summary:\n")
	print(data.frame(rbind(unlist(spweights.constants(x,
		zero.policy=zero.policy))[c(1, 5:8)]), row.names=style))
	invisible(x)

}

