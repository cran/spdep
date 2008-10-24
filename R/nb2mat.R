# Copyright 2001-7 by Roger Bivand, Markus Reder and Werner Mueller
#


nb2mat <- function(neighbours, glist=NULL, style="W", zero.policy=FALSE)
{
	if(!inherits(neighbours, "nb")) stop("Not a neighbours list")
	listw <- nb2listw(neighbours, glist=glist, style=style,
		zero.policy=zero.policy)
	res <- listw2mat(listw)
	attr(res, "call") <- match.call()
	res
}

listw2mat <- function(listw) {
	n <- length(listw$neighbours)
	if (n < 1) stop("non-positive number of entities")
	cardnb <- card(listw$neighbours)
	if (any(is.na(unlist(listw$weights))))
		stop ("NAs in general weights list")
	res <- matrix(0, nrow=n, ncol=n)
	for (i in 1:n)
	    if (cardnb[i] > 0)
		res[i, listw$neighbours[[i]]] <- listw$weights[[i]]
	if (!is.null(attr(listw, "region.id")))
		row.names(res) <- attr(listw, "region.id")
	res
}

invIrM <- function(neighbours, rho, glist=NULL, style="W", method="solve", 
	feasible=NULL) {
	if(class(neighbours) != "nb") stop("Not a neighbours list")
	invIrW(nb2listw(neighbours, glist=glist, style=style), rho=rho, 
		method=method, feasible=feasible)
}

invIrW <- function(listw, rho, method="solve", feasible=NULL) {
	if(!inherits(listw, "listw")) stop("Not a weights list")
	n <- length(listw$neighbours)
	V <- listw2mat(listw)
	if (is.null(feasible) || (is.logical(feasible) && !feasible)) {
		V <- listw2mat(listw)
		e <- eigen(V, only.values = TRUE)$values
		if (is.complex(e)) feasible <- 1/(range(Re(e)))
		else feasible <- 1/(range(e))
		if (rho <= feasible[1] || rho >= feasible[2])
			stop(paste("Rho outside feasible range:", feasible))
	}
	if (method == "chol"){
		if (listw$style %in% c("W", "S") && !(can.be.simmed(listw)))
			stop("Cholesky method requires symmetric weights")
		if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
			stop("Cholesky method requires symmetric weights")
		if (listw$style %in% c("W", "S")) {
			V <- listw2mat(listw2U(similar.listw(listw)))
		}
		mat <- diag(n) - rho * V
		res <- chol2inv(chol(mat))
	} else if (method == "solve") {
		mat <- diag(n) - rho * V
		res <- solve(mat)
	} else stop("unknown method")
	attr(res, "call") <- match.call()
	res
}

mat2listw <- function(x, row.names=NULL) {
	if (!is.matrix(x)) stop("x is not a matrix")
	n <- nrow(x)
	if (n < 1) stop("non-positive number of entities")
	m <- ncol(x)
	if (n != m) stop("x must be a square matrix")
	if (any(x < 0)) stop("values in x cannot be negative")
	if (any(is.na(x))) stop("NA values in x not allowed")
    	if (!is.null(row.names)) {
		if(length(row.names) != n)
            		stop("row.names wrong length")
		if (length(unique(row.names)) != length(row.names))
	    		stop("non-unique row.names given")
    	}
    	if (is.null(row.names)) {
		if (!is.null(row.names(x))) {
			row.names <- row.names(x)
		} else {
			row.names <- as.character(1:n)
		}
	}
	style <- "M"
	neighbours <- vector(mode="list", length=n)
	weights <- vector(mode="list", length=n)
	for (i in 1:n) {
		nbs  <- which(x[i,] > 0.0)
		if (length(nbs) > 0) {
			neighbours[[i]] <- nbs
			weights[[i]] <- as.double(x[i, nbs]) # Laurajean Lewis
		} else {
			neighbours[[i]] <- as.integer(0)
		}
	}
	attr(weights, "mode") <- "unknown" # Brian Rubineau
	class(neighbours) <- "nb"
	attr(neighbours, "region.id") <- row.names
 	attr(neighbours, "call") <- NA
        attr(neighbours, "sym") <- is.symmetric.nb(neighbours, 
		verbose=FALSE, force=TRUE)
	res <- list(style=style, neighbours=neighbours, weights=weights)
	class(res) <- c("listw", "nb")
	attr(res, "region.id") <- attr(neighbours, "region.id")
	attr(res, "call") <- match.call()
	res
}
