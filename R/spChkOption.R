# Copyright 2003 by Roger Bivand 

set.spChkOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("spChkID", env = .spChkOption)
	assign("spChkID", check, env = .spChkOption)
	res
}

get.spChkOption <- function() {
	get("spChkID", env = .spChkOption)
}

chkIDs <- function(x, listw) {
	if (!is.array(x) & !is.data.frame(x)) {
	    if (is.null(xn <- names(x)))
	    	stop(paste(deparse(substitute(x)), "has no names"))
	} else {
	    if (is.null(xn <- rownames(x)))
		stop(paste(deparse(substitute(x)), "has no row names"))
	}
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (is.null(ln <- attr(listw, "region.id")))
		stop(paste(deparse(substitute(listw)), "has no region IDs"))
	if (length(ln) != length(xn)) stop("objects of different length")
	res <- all(ln == xn)
	res
}

spNamedVec <- function(var, data) {
	if (!is.array(data) & !is.data.frame(data))
		stop(paste(deparse(substitute(data)),
			"not an array or data frame"))
	if (!is.character(var) & !is.numeric(var)) 
		stop("variable name wrong type") 
	res <- try(data[,var])
	if (inherits(res, "try-error")) 
		stop(paste(deparse(substitute(var)), "not found"))
	nms <- rownames(data)
	if (is.null(nms)) nms <- as.character(1:length(res))
	names(res) <- nms
	res
}
