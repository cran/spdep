# Copyright 2001-3 by Roger Bivand 
#

spweights.constants <- function(listw, zero.policy=FALSE) {
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	cards <- card(listw$neighbours)
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
	n <- length(which(cards > 0))
	n1 <- n - 1
	n2 <- n - 2
	n3 <- n - 3
	nn <- n*n
	S0 <- Szero(listw)
	S1 <- 0
	S2 <- 0
	for (i in 1:length(listw$neighbours)) {
		cond <- TRUE
		if (zero.policy && cards[i] == 0) cond <- FALSE
		if (cond) {
			if (cards[i] == 0)
				stop(paste("region", i,
					"has no neighbours"))
			ij <- listw$neighbours[[i]]
			wij <- listw$weights[[i]]
			dm0 <- 0
			dm1 <- 0
			for (j in 1:length(ij)) {
				dij <- wij[j]
				ij.j <- ij[j]
				ij.lkup <- which(listw$neighbours[[ij.j]] == i)
				if (length(ij.lkup) == 1)
					dji <- listw$weights[[ij.j]][ij.lkup]
				else dji <- 0
				dm0 <- dm0 + dij
				dm1 <- dm1 + dji
				S1 <- S1 + (dij + dji)^2
			}
			S2 <- S2 + (dm0 + dm1)^2
		}
	}
	S1 <- S1 * 0.5
	invisible(list(n=n, n1=n1, n2=n2, n3=n3, nn=nn, S0=S0, S1=S1, S2=S2))
}

Szero <- function(listw) {
	sum(unlist(listw$weights))
}

lag.listw <- function(listw, x, zero.policy=FALSE) {
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(listw)),
		"not numeric"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	if (is.vector(x)) {
		if (length(x) != n) stop("object lengths differ")
		res <- .Call("lagw", listw$neighbours, listw$weights,
			as.double(x), as.integer(cardnb),
			as.logical(zero.policy), PACKAGE="spdep")
	} else if (is.matrix(x)) {
		if (nrow(x) != n) stop("object lengths differ")
		res <- matrix(0, nrow=nrow(x), ncol=ncol(x))
		for (i in 1:ncol(x)) {
			res[,i] <- .Call("lagw", listw$neighbours,
				listw$weights, as.double(x[,i]),
				as.integer(cardnb), as.logical(zero.policy),
				PACKAGE="spdep")

		}
	} else {
		stop(paste(deparse(substitute(x)),
			"neither a numeric vector or matrix"))
	}
	if (any(is.na(res))) warning("NAs in lagged values")
	invisible(res)
}

listw2U <- function(listw) {
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	nb <- listw$neighbours
	wts <- listw$weights
	style <- paste(listw$style, "U", sep="")
	sym <- is.symmetric.nb(nb, FALSE, TRUE)
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	nlist <- vector(mode="list", length=n)
	attr(nlist, "region.id") <- attr(nb, "region.id")
	class(nlist) <- "nb"
	vlist <- vector(mode="list", length=n)
	attr(vlist, as.character(style)) <- TRUE
	if (sym) {
		nlist <- vector(mode="list", length=n)
		attr(nlist, "region.id") <- attr(nb, "region.id")
		class(nlist) <- "nb"
		for (i in 1:n) {
			inb <- nb[[i]]
			nlist[[i]] <- inb
			iwt <- wts[[i]]
			icd <- cardnb[i]
			if (icd > 0) {
			    for (j in 1:icd) {
				vlist[[i]][j] <- 0.5 *
				(iwt[j]+wts[[inb[j]]][which(nb[[inb[j]]] == i)])
			    }
			}
		}
	} else {
		nlist <- make.sym.nb(nb)
		for (i in 1:n) {
			inb <- nb[[i]]
			inl <- nlist[[i]]
			if (inl > 0) {
			    iwt <- wts[[i]]
			    vlist[[i]] <- numeric(length=length(inl))
			    for (j in 1:length(inl)) {
				if (inl[j] %in% inb) a <-
					iwt[which(inb == inl[j])]
				else a <- 0
				if (i %in% nb[[inl[j]]]) b <-
					wts[[inl[j]]][which(nb[[inl[j]]] == i)]
				else b <- 0
				vlist[[i]][j] <- 0.5 * (a + b)
			    }
			}
		}
	}
	res <- list(style=style, neighbours=nlist, weights=vlist)
	class(res) <- "listw"
	attr(res, "region.id") <- attr(nb, "region.id")
	attr(res, "call") <- match.call()
	attr(res, "U") <- TRUE
	invisible(res)
}


listw2star <- function(listw, ireg, style, n, D, a, zero.policy) {
    nb <- vector(mode="list", length=n)
    class(nb) <- "nb"
    wts <- vector(mode="list", length=n)
    for (i in 1:n) nb[[i]] <- 0
    inb <- listw$neighbours[[ireg]]
    iwts <- listw$weights[[ireg]]
    cond <- TRUE
    if (inb == 0 || length(inb) == 0 || is.null(iwts)) cond <- FALSE
    if (!cond && !zero.policy) stop("No-neighbour region found")
    if (style == "W") iwts <- (n*D[ireg]*iwts) / 2
    else if (style == "S") iwts <- ((n^2)*D[ireg]*iwts) / (2*a)
    else if (style == "C") iwts <- ((n^2)*iwts) / (2*a)
    if (cond) {
    	nb[[ireg]] <- inb
    	wts[[ireg]] <- iwts
    	for (j in 1:length(inb)) {
            jj <- inb[j]
            nb[[jj]] <- ireg
            wts[[jj]] <- iwts[j]
	}
    }
    res <- list(style=style, neighbours=nb, weights=wts)
    class(res) <- c("listw", "star")
    attr(res, "region.id") <- attr(listw, "region.id")
    res
}


