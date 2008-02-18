# Copyright 2001-2006 by Roger Bivand 
#
	


poly2nb <- function(pl, row.names=NULL, snap=sqrt(.Machine$double.eps),
	queen=TRUE) {
	if (!inherits(pl, "polylist")) {
		if (extends(class(pl), "SpatialPolygons"))
			pl <- maptools:::.SpP2polylist(pl)
		else stop("Not a polygon list")
	}
	if (inherits(pl, "multiparts")) stop("Convert to newer polylist format")
	n <- length(pl)
	if (n < 1) stop("non-positive number of entities")
	if (is.null(row.names)) regid <- attr(pl, "region.id")
	else regid <- NULL
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
		if (n < 1) stop("non-positive number of entities")
		res <- matrix(0, nrow=n, ncol=4)
		for (i in 1:n) res[i,] <- attr(pl[[i]], "bbox")
		res
	}
	bb <- poly2bbs(pl)
	if (storage.mode(bb) != "double") storage.mode(bb) <- "double"
	dsnap <- as.double(snap)
	bb[,1] <- bb[,1] - dsnap
	bb[,2] <- bb[,2] - dsnap
	bb[,3] <- bb[,3] + dsnap
	bb[,4] <- bb[,4] + dsnap
	nrs <- integer(n)
	for (i in 1:n) {
		pl[[i]] <- na.omit(pl[[i]][-1,])
		nrs[i] <- as.integer(nrow(pl[[i]]))
		pl[[i]] <- as.double(pl[[i]])
	}
	


	polypoly2 <- function(poly1, nrs1, poly2, nrs2, snap) {
		if (any(nrs1 == 0 || nrs2 == 0)) return(as.integer(0))
		res <- .Call("polypoly", poly1, nrs1, poly2, 
			nrs2, snap, PACKAGE="spdep")
		res
	}

	ans <- vector(mode="list", length=n)
	for (i in 1:n) ans[[i]] <- integer(0)
	criterion <- ifelse(queen, 0, 1)
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			jhit <- .Call("spInsiders", bb[i,], 
				bb[j,], PACKAGE="spdep")
			if (jhit > 0) {
			    khit <- 0
			    khit <- polypoly2(pl[[i]], nrs[i], pl[[j]], 
				nrs[j],dsnap)

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
	ans <- sym.attr.nb(ans)
	ans
}


