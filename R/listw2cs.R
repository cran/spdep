# Copyright 2004 Roger Bivand


asMatrixCsrListw <- function(listw, zero.policy=FALSE) {
	if(!inherits(listw, "listw")) stop("not a listw object")
	n <- length(listw$neighbours)
	ra <- unlist(listw$weights)
	ja <- unlist(listw$neighbours)
# omitted zero policy found at Kohren Sahlis 2005-11-24
        if (zero.policy) {
		zeros <- which(ja == 0)
		if (length(zeros) > 0) ja <- ja[-zeros]
	}
        if (length(ja) != length(ra))
	  stop("different numbers of weights and numbers - wrong zero policy?")
	ia <- as.integer(cumsum(c(1, card(listw$neighbours))))
	res <- new("matrix.csr", ra=ra, ja=ja, ia=ia, dimension=c(n,n))
	res
}

asMatrixCsrI <- function(n) {
	if (n < 1) stop("matrix must have positive dimensions")
	I <- new("matrix.csr", ra=rep(1, n), ja=1:n, ia=1:(n+1), 
		dimension=c(n,n))
	I
}
asMatrixCsrIrW <- function(W, rho) {
	if(!inherits(W, "matrix.csr")) stop("not a matrix.csr object")
	n <- dim(W)[1]
	I <- asMatrixCsrI(n)
	IrW <- I - rho * W
	IrW
}

asListwMatrixCsr <- function(mcsr) {
	if(!is.matrix.csr(mcsr)) stop("not a matrix.csr object")
	dim <- mcsr@dimension
	if (dim[1] != dim[2]) warning("rectangular matrix")
	n <- dim[1]
	if (n < 1) stop("non-positive dimension")
	ra <- mcsr@ra
	ja <- mcsr@ja
	ia <- mcsr@ia
	dia <- diff(mcsr@ia)
	if(length(dia) != n) stop("dimension does not match row indices")
	if (any(dia < 1) | any(dia > n)) stop("row indices out of range")
	if (any(ja < 1) | any(ja > n)) stop("column indices out of range")
	region.id <- as.character(1:n)
	nlist <- vector(mode="list", length=n)
	class(nlist) <- "nb"
	attr(nlist, "region.id") <- region.id
	vlist <- vector(mode="list", length=n)
	for (i in 1:n) {
		if (dia[i] > 0) {
			interval <- ia[i]:(ia[i]+dia[i]-1)
			nlist[[i]] <- as.integer(ja[interval])
			vlist[[i]] <- as.double(ra[interval])
		} else {
			nlist[[i]] <- as.integer(0)
		}
	}
	res <- list(style=as.character(NA), neighbours=nlist, weights=vlist)
	class(res) <- c("listw", "nb")
	attr(res, "region.id") <- region.id
	attr(res, "call") <- match.call()
	res
}


