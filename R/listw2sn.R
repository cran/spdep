# Copyright 2001 by Roger Bivand
#

listw2sn <- function(listw) {
	if(!inherits(listw, "listw")) stop("not a listw object")
	z <- .Call("listw2sn", listw$neighbours, listw$weights,
		PACKAGE="spdep")
	res <- as.data.frame(list(from=z[[1]], to=z[[2]], weights=z[[3]]))
	class(res) <- c(class(res), "spatial.neighbour")
	attr(res, "region.id") <- attr(listw, "region.id")
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(res, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(res, "weights.attrs") <- weights.attrs
	attr(res, "listw.call") <- attr(listw, "call")
	invisible(res)
}
