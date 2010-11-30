# Copyright 2003-2010 by Roger Bivand 

set.spChkOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("spChkID", env = .spdepOptions)
	assign("spChkID", check, env = .spdepOptions)
	res
}

get.spChkOption <- function() {
	get("spChkID", env = .spdepOptions)
}

set.VerboseOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("verbose", env = .spdepOptions)
	assign("verbose", check, env = .spdepOptions)
	res
}

get.VerboseOption <- function() {
	get("verbose", env = .spdepOptions)
}

set.ZeroPolicyOption <- function(check) {
	if (!is.logical(check)) stop ("logical argument required")
	res <- get("zeroPolicy", env = .spdepOptions)
	assign("zeroPolicy", check, env = .spdepOptions)
	res
}

get.ZeroPolicyOption <- function() {
	get("zeroPolicy", env = .spdepOptions)
}

set.ClusterOption <- function(cl) {
	if (!is.null(cl)) {
            if (!inherits(cl, "cluster")) 
                stop ("cluster required")
            clusterEvalQ(cl, library(spdep))
        }
        if (is.null(cl)) clusterEvalQ(get.ClusterOption(),
            detach(package:spdep))
	assign("cl", cl, env = .spdepOptions)
        invisible(NULL)
}

get.ClusterOption  <- function() {
	get("cl", env = .spdepOptions)
}

get.rlecuyerSeedOption  <- function() {
	get("rlecuyerSeed", env = .spdepOptions)
}

set.rlecuyerSeedOption  <- function(seed) {
    if (length(seed) != 6) stop("Six integer values required")
    if (storage.mode(seed) != "integer") seed <- as.integer(seed)
    assign("rlecuyerSeed", seed, env = .spdepOptions)
    invisible(NULL)
}

chkIDs <- function (x, listw) 
{
    if (!is.array(x) & !is.data.frame(x)) {
        if (is.null(xn <- names(x))) 
            stop(paste(deparse(substitute(x)), "has no names"))
    }
    else {
        if (is.null(xn <- rownames(x))) 
            stop(paste(deparse(substitute(x)), "has no row names"))
    }
    if (!inherits(listw, "nb")) 
        stop(paste(deparse(substitute(listw)), "is not an listw  or nb object"))
    if (is.null(ln <- attr(listw, "region.id"))) 
        stop(paste(deparse(substitute(listw)), "has no region IDs"))
    if (length(ln) != length(xn)) 
        stop("objects of different length")
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
