# Copyright 2001 by Roger Bivand
#


nb2mat <- function(neighbours, glist=NULL, style="W", zero.policy=FALSE)
{
	if(class(neighbours) != "nb") stop("Not a neighbours list")
	listw <- nb2listw(neighbours, glist=glist, style=style,
		zero.policy=zero.policy)
	res <- listw2mat(listw)
	attr(res, "call") <- match.call()
	invisible(res)
}

listw2mat <- function(listw) {
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	if (any(is.na(unlist(listw$weights))))
		stop ("NAs in general weights list")
	res <- matrix(0, nrow=n, ncol=n)
	for (i in 1:n)
	    if (cardnb[i] > 0)
		res[i, listw$neighbours[[i]]] <- listw$weights[[i]]
	invisible(res)
}

invIrM <- function(neighbours, rho, glist=NULL, style="W") {
	if(class(neighbours) != "nb") stop("Not a neighbours list")
	n <- length(neighbours)
	V <- nb2mat(neighbours, glist, style)
	feasible <- 1/(range(eigen(V, only.values=TRUE)$values))
	if (rho <= feasible[1] || rho >= feasible[2])
		stop(paste("Rho outside feasible range:", feasible))
	mat <- diag(n) - rho * V
	res <- solve(mat)
	attr(res, "call") <- match.call()
	invisible(res)
}
