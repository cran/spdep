# Copyright 2001 by Roger S. Bivand. 
#

nb2listw <- function(neighbours, glist=NULL, style="W", zero.policy=FALSE)
{
	if(class(neighbours) != "nb") stop("Not a neighbours list")
	n <- length(neighbours)
	cardnb <- card(neighbours)
	if (!zero.policy)
		if (any(cardnb == 0)) stop("Empty neighbour sets found")
	vlist <- vector(mode="list", length=n)
	if (is.null(glist)) {
		glist <- vector(mode="list", length=n)
		for (i in 1:n)
			if(cardnb[i] > 0) glist[[i]] <- rep(1, length=cardnb[i])
		attr(vlist, "binary") <- TRUE
	} else {
		attr(vlist, "general") <- TRUE
		source <- deparse(substitute(glist))
		attr(vlist, as.character(source)) <- TRUE
	}
	if (length(glist) != n) stop("glist wrong length")
	if (cardnb != unlist(lapply(glist, length)))
		stop("neighbours and glist do not conform")
	if (any(is.na(unlist(glist))))
		stop ("NAs in general weights list")
	attr(vlist, as.character(style)) <- TRUE
	if (style == "W") {
		d <- unlist(lapply(glist, sum))
		for (i in 1:n) {
			if (cardnb[i] > 0) vlist[[i]] <- (1/d[i]) * glist[[i]]
		}
	}
	if (style == "B") {
		for (i in 1:n) {
			if (cardnb[i] > 0) vlist[[i]] <- rep(1, cardnb[i])
		}
	}
	if (style == "C") {
		D <- sum(unlist(glist))
		if (is.na(D) || !(D > 0))
			stop(paste("Failure in sum of weights:", D))
		if (zero.policy) {
			eff.n <- n - length(which(cardnb == 0))
			if (eff.n < 1) stop("No valid observations")
		} else eff.n <- n
		for (i in 1:n) {
			if (cardnb[i] > 0)
				vlist[[i]] <- (eff.n/D) * glist[[i]]
		}
	}
	if (style == "S") {
		glist2 <- lapply(glist, function(x) x^2)
		q <- sqrt(unlist(lapply(glist2, sum)))
		for (i in 1:n) {
			if (cardnb[i] > 0)
				glist[[i]] <- (1/q[i]) * glist[[i]]
		}
		Q <- sum(unlist(glist))
		if (is.na(Q) || !(Q > 0))
		    stop(paste("Failure in sum of intermediate weights:", Q))
		if (zero.policy) {
			eff.n <- n - length(which(cardnb == 0))
			if (eff.n < 1) stop("No valid observations")
		} else eff.n <- n
		for (i in 1:n) {
			if (cardnb[i] > 0)
				vlist[[i]] <- (eff.n/Q) * glist[[i]]
		}
	}
	style <- style
	if (!zero.policy)
		if (any(is.na(unlist(vlist))))
			stop ("NAs in coding scheme weights list")
	res <- list(style=style, neighbours=neighbours, weights=vlist)
	class(res) <- "listw"
	attr(res, "region.id") <- attr(neighbours, "region.id")
	attr(res, "call") <- match.call()
	invisible(res)
}

