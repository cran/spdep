# Copyright 2001-2004 by Roger Bivand 
#
	


poly2nb <- function(pl, row.names=NULL, snap=sqrt(.Machine$double.eps),
	queen=TRUE) {
	if (!inherits(pl, "polylist")) stop("Not a polygon list")
	if (inherits(pl, "multiparts")) stop("Convert to newer polylist format")
	n <- length(pl)
	regid <- attr(pl, "region.id")
	if (is.null(regid)) {
		if(is.null(row.names)) regid <- as.character(1:n)
		else {
			if(length(row.names) != n)
				stop("row.names wrong length")
			else if (length(unique(row.names)) != length(row.names))
	    			stop("non-unique row.names given")
			else regid <- row.names
		}
	}
	poly2bbs <- function(pl) {
		n <- length(pl)
		res <- matrix(0, nrow=n, ncol=4)
		for (i in 1:n) res[i,] <- attr(pl[[i]], "bbox")
		res
	}
	bb <- poly2bbs(pl)
	


	polypoly2 <- function(poly1, poly2, snap) {
		n1 <- nrow(poly1)
		n2 <- nrow(poly2)
		res <- .Call("polypoly", as.double(poly1), 
			as.integer(n1), as.double(poly2), 
			as.integer(n2), as.double(snap), PACKAGE="spdep")
		res
	}

	ans <- vector(mode="list", length=n)
	for (i in 1:n) ans[[i]] <- integer(0)
	criterion <- ifelse(queen, 0, 1)
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			jhit <- .Call("spInsiders", as.double(bb[i,]), 
				as.double(bb[j,]), PACKAGE="spdep")
			if (jhit > 0) {
			    khit <- 0
			    khit <- polypoly2(na.omit(pl[[i]]), 
				na.omit(pl[[j]]), snap)

			    if (khit > criterion) {
				ans[[i]] <- c(ans[[i]], j)
				ans[[j]] <- c(ans[[j]], i)
			    }
			}
		}
	}
	for (i in 1:n) ans[[i]] <- sort(ans[[i]])
	class(ans) <- "nb"
	attr(ans, "region.id") <- regid
	attr(ans, "call") <- match.call()
	if (queen) attr(ans, "type") <- "queen"
	else attr(ans, "type") <- "rook"
	invisible(ans)
}


