# Copyright 2001-2003 by Roger Bivand
#


plotpolys <- function(pl, bb, col=NA, border=par("fg"), add=FALSE, 
	xlim=NULL, ylim=NULL) {
	if (all(class(pl) != "polylist")) stop("Not a polygon list")
	if (!add) {
		if (is.null(xlim)) xlim <- c(min(bb[,1]), max(bb[,3]))
		if (is.null(ylim)) ylim <- c(min(bb[,2]), max(bb[,4]))
		plot(x=bb[,1], y=bb[,4], xlim=xlim, ylim=ylim, type="n",
		asp=1, xlab="", ylab="")
	}
	if (length(col) != length(pl)) {
		col <- rep(col, length(pl), length(pl))
	}
	for (j in 1:length(pl)) {
		if ("multiparts" %in% class(pl)) {
			for (k in 1:length(pl[[j]]))
				polygon(pl[[j]][[k]], col=col[j], border=border)
		} else {
			polygon(pl[[j]], col=col[j], border=border)
		}
	}
}
	
poly2nb <- function(pl, bb, row.names=NULL, snap=sqrt(.Machine$double.eps),
	queen=TRUE) {
	if (all(class(pl) != "polylist")) stop("Not a polygon list")
	if ("multiparts" %in% class(pl)) stop("No multiparts yet")
	n <- length(pl)
	regid <- attr(pl, "region.id")
	if (is.null(regid)) {
		if(is.null(row.names)) regid <- as.character(1:n)
		else {
			if(length(row.names) != n)
				stop("row.names wrong length")
			else regid <- row.names
		}
	}
	if (nrow(bb) != n)
		stop("Number of polygons not equal to number of bounding boxes")

	between <- function(x, low, up) {return(x >= low && x <= up)}

	pipbb <- function(pt, bb) {
		return(between(pt[1], bb[1], bb[3]) && 
			between(pt[2], bb[2], bb[4]))
	}

	polypoly <- function(poly1, poly2, snap) {
		snap2 <- snap^2
		a <- outer(poly1[,1], poly2[,1], "-")
		b <- outer(poly1[,2], poly2[,2], "-")
		c <- a^2 + b^2
		d <- rle(unlist(apply(c,2,function(x) {
			length(which(x < snap2)) > 0} )))
		res <- length(which(d$values == TRUE))
		if (res > 0) res <- sum(d$length[d$values == TRUE])
		res
	}

	pipbbij <- function(bbi, bbj) {
		jhit <- logical(8)
	    	jhit[1] <- pipbb(c(bbi[1], bbi[2]), bb[j,])
		jhit[2] <- pipbb(c(bbi[1], bbi[4]), bb[j,])
		jhit[3] <- pipbb(c(bbi[3], bbi[2]), bb[j,])
		jhit[4] <- pipbb(c(bbi[3], bbi[4]), bb[j,])
		jhit[5] <- pipbb(c(bb[j,1], bb[j,2]), bbi)
		jhit[6] <- pipbb(c(bb[j,1], bb[j,4]), bbi)
		jhit[7] <- pipbb(c(bb[j,3], bb[j,2]), bbi)
		jhit[8] <- pipbb(c(bb[j,3], bb[j,4]), bbi)
		return(jhit)
	}

	ans <- vector(mode="list", length=n)
	for (i in 1:n) ans[[i]] <- integer(0)
	criterion <- ifelse(queen, 0, 1)
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			jhit <- pipbbij(bb[i,], bb[j,])
			if (any(jhit)) {
			    khit <- 0
			    khit <- polypoly(pl[[i]], pl[[j]], snap)
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

