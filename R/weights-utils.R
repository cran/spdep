# Copyright 2001-4 by Roger Bivand 
#


is.symmetric.nb <- function(nb, verbose=TRUE, force=FALSE)
{
	if(!inherits(nb, "nb")) stop("Not neighbours list")
	nbsym <- attr(nb, "sym")
	if(!is.null(nbsym)) res <- nbsym
	if(force || is.null(nbsym)) {
		res <- .Call("symtest", nb=nb, card=as.integer(card(nb)),
			verbose=as.logical(verbose), PACKAGE="spdep")
	}
	if(!res && verbose) cat("Non-symmetric neighbours list\n")
	invisible(res)
}

sym.attr.nb <- function(nb) {
	if(!inherits(nb, "nb")) stop("Not neighbours list")
	nbsym <- attr(nb, "sym")
	if(is.null(nbsym))
		attr(nb, "sym") <- is.symmetric.nb(nb, verbose=FALSE,
			force=TRUE)
	invisible(nb)
}

include.self <- function(nb) {
	if (!is.null(attributes(nb)$self.included) &&
		(as.logical(attributes(nb)$self.included)))
		stop("Self already included")
	n <- length(nb)
	nc <- card(nb)
	for (i in 1:n) {
		if (nc[i] > 0) {
			nb[[i]] <- sort(c(i, nb[[i]]))
		} else {
			nb[[i]] <- i
		}
	}
		
	attr(nb, "self.included") <- TRUE
	invisible(nb)
}

# Copyright 2001 by Nicholas Lewin-Koh 

make.sym.nb <- function(nb){
	if(!inherits(nb, "nb")) stop("Not neighbours list")
	if (is.symmetric.nb(nb, FALSE, TRUE)) {
		res <- nb
	} else {
        	k <- unlist(lapply(nb,length))
        	to <- unlist(nb)
        	from <- NULL
        	res <- vector(mode="list", length=length(nb))
        	for(i in 1:length(nb)){
        		from <- c(from,rep(i,k[i]))
        	}
        	for(i in 1:length(nb)){
        		res[[i]] <- sort(unique(c(to[from==i],from[to==i])))
        		if(length(res[[i]]) == 0) res[[i]] <- 0
        	}
        	attr(res, "region.id") <- attr(nb,"region.id")
        	attr(res, "call") <- attr(nb, "call")
        	attr(res, "type") <- attr(nb, "type")
        	attr(res, "sym") <- TRUE
        	class(res) <- "nb"
	}
	invisible(res)
}

