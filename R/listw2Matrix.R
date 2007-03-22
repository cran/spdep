# Copyright 2006-7 by Roger Bivand
#

as_dgRMatrix_listw <- function(listw) {
	if(!inherits(listw, "listw")) stop("not a listw object")
	n <- length(listw$neighbours)
	cardw <- card(listw$neighbours)
	p0 <- as.integer(c(0, cumsum(cardw)))
	scard <- sum(cardw)
	z <- .Call("listw2dgR", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard), PACKAGE="spdep")
	res <- new("dgRMatrix", j=z[[1]], p=p0, Dim=as.integer(c(n, n)),
		x=z[[2]])
	res
}

as_dsTMatrix_listw <- function(listw) {
	if (!inherits(listw, "listw")) stop("not a listw object")
	if (!is.symmetric.glist(listw$neighbours, listw$weights))
		stop("not a symmetric matrix")
	n <- length(listw$neighbours)
	cardw <- card(listw$neighbours)
	scard <- sum(cardw)
	if (scard %% 2 != 0) stop("odd neighbours sum")
	z <- .Call("listw2dsT", listw$neighbours, listw$weights,
		as.integer(cardw), as.integer(scard/2), PACKAGE="spdep")

	res <- new("dsTMatrix", i=z[[1]], j=z[[2]], Dim=as.integer(c(n, n)),
		x=z[[3]])
	res
}

as_dgCMatrix_I <- function(n) {
	if (n < 1) stop("matrix must have positive dimensions")
	I <- as(Diagonal(n), "sparseMatrix")
	I
}

as_dgCMatrix_IrW <- function(W, rho) {
	if(!inherits(W, "dsTMatrix")) stop("not a dsTMatrix object")
	n <- dim(W)[1]
	I <- as_dgCMatrix_I(n)
	IrW <- I - rho * W
	IrW
}

Jacobian_W <- function(W, rho) {
	IrW <- as_dgCMatrix_IrW(W, rho)
	logdet <- sum(2*log(diag(chol(as(IrW, "dsCMatrix")))))
	logdet
}

