# Copyright 2001 by Roger Bivand
#

subset.nb <- function(x, subset, ...) {
    if (class(x) != "nb") stop("not a neighbours list")
    if (!is.logical(subset)) stop("subset not a logical vector")
    n <- length(x)
    if (n != length(subset))
	stop("neighours list and subset vector different lengths")
    old.ids <- 1:n
    new.ids <- match(old.ids, which(subset))
    reg.id <- subset.default(attr(x, "region.id"), subset)
    x <- sym.attr.nb(x)
    xattrs <- names(attributes(x))
    z <- subset.default(x, subset)
    nz <- length(z)
    for (i in 1:nz) {
	zi <- z[[i]]
	res <- NULL
	for (j in 1:length(zi)) {
	    a <- new.ids[zi[j]]
	    if (!is.na(a)) res <- c(res, a)
	}
	if (is.null(res)) z[[i]] <- 0
	else z[[i]] <- sort(unique(res))
    }
    attr(z, "region.id") <- reg.id
    for (i in 1:length(xattrs)) {
	if (xattrs[i] != "region.id")
	    attr(z, xattrs[i]) <- attr(x, xattrs[i])
    }
    invisible(z)
}

