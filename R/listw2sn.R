# Copyright 2001-3 by Luc Anselin and Roger Bivand
#

listw2sn <- function(listw) {
	if(!inherits(listw, "listw")) stop("not a listw object")
	z <- .Call("listw2sn", listw$neighbours, listw$weights,
		PACKAGE="spdep")
	res <- as.data.frame(list(from=z[[1]], to=z[[2]], weights=z[[3]]))
	class(res) <- c(class(res), "spatial.neighbour")
	attr(res, "n") <- length(attr(listw, "region.id"))
	attr(res, "region.id") <- attr(listw, "region.id")
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(res, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(res, "weights.attrs") <- weights.attrs
	if (!(is.null(attr(listw, "GeoDa"))))
		attr(res, "GeoDa") <- attr(listw, "GeoDa")
	attr(res, "listw.call") <- attr(listw, "call")
	invisible(res)
}

sn2listw <- function(sn) {
	if(!inherits(sn, "spatial.neighbour")) 
	    stop("not a spatial.neighbour object")
	n <- attr(sn, "n")
	region.id <- attr(sn, "region.id")
	nlist <- vector(mode="list", length=n)
	vlist <- vector(mode="list", length=n)
	rle.sn <- rle(sn[,1])
	if (n != length(rle.sn$lengths))
	    warning(paste(paste(region.id[which(!(1:n %in% rle.sn$values))], 
		collapse=", "), "are not origins"))
	cs1.sn <- cumsum(rle.sn$lengths)
	cs0.sn <- c(1, cs1.sn[1:(n-1)]+1)
	ii <- 1
	for (i in 1:n) {
		if (rle.sn$value[ii] == i) {
			nlist[[i]] <- as.integer(sn[cs0.sn[ii]:cs1.sn[ii],2])
			vlist[[i]] <- as.double(sn[cs0.sn[ii]:cs1.sn[ii],3])
			ii <- ii+1
		} else {
			nlist[[i]] <- as.integer(0)
		}
	}
	res <- list(style=as.character(NA), neighbours=nlist, weights=vlist)
	class(res) <- c("listw", "nb")
	if (!(is.null(attr(sn, "GeoDa"))))
		attr(res, "GeoDa") <- attr(sn, "GeoDa")
	attr(res, "region.id") <- region.id
	attr(res, "call") <- match.call()
	invisible(res)
}

# LA 6/28/03 read.gwt
# LA 7/12/03 revised, sorted ids

read.gwt2nb <- function(file, region.id=NULL) {
	con <- file(file, open="r")   #opens the file
	firstline <- unlist(strsplit(readLines(con,1)," "))
	if (length(firstline) == 4) {
		n <- as.integer(firstline[2])
		shpfile <- firstline[3]
		ind <- firstline[4]
		if (ind != deparse(substitute(region.id)))
			warning(paste("region.id not named", ind))
	} else if (length(firstline) == 1) {
		n <- as.integer(firstline[1])
		shpfile <- as.character(NA)
		ind <- as.character(NA)
		warning("Old-style GWT file")
	} else stop("Invalid header line format for GWT file")
	close(con)
	nseq <- 1:n
	if (is.null(region.id)) region.id <- nseq
	if (n != length(region.id))
		stop("Mismatch in dimensions of GWT file and region.id")
	odij <- read.table(file, skip=1)
	qorder <- order(odij[,1],odij[,2])
	odij <- odij[qorder,]
	origvec <- unique(odij[,1])
	if (!all(nseq %in% origvec))
		warning(paste(paste(region.id[which(!(nseq %in% origvec))], 
			collapse=", "), "are not origins"))
	destvec <- unique(odij[,2])
	if (!all(nseq %in% destvec))
		warning(paste(paste(region.id[which(!(nseq %in% destvec))], 
			collapse=", "), "are not destinations"))

	res <- vector(mode="list", length=n)
	vlist <- vector(mode="list", length=n)
	rle.sn <- rle(odij[,1])
	cs1.sn <- cumsum(rle.sn$lengths)
	cs0.sn <- c(1, cs1.sn[1:(n-1)]+1)
	ii <- 1
	for (i in 1:n) {
		if (rle.sn$value[ii] == i) {
			res[[i]] <- as.integer(odij[cs0.sn[ii]:cs1.sn[ii],2])
			vlist[[i]] <- as.double(odij[cs0.sn[ii]:cs1.sn[ii],3])
			ii <- ii+1
		} else {
			res[[i]] <- as.integer(0)
		}
	}

	class(res) <- c("nb", "GWT")
	attr(res, "region.id") <- region.id
	attr(res, "neighbours.attrs") <- as.character(NA)
	attr(res, "weights.attrs") <- as.character(NA)
	attr(res, "GeoDa") <- list(dist=vlist, shpfile=shpfile, ind=ind)
	attr(res, "call") <- match.call()
	attr(res, "n") <- n
	invisible(res)
}

write.sn2gwt <- function(sn, file, shpfile=NULL, ind=NULL) {
	if(!inherits(sn, "spatial.neighbour")) 
	    stop("not a spatial.neighbour object")
	n <- attr(sn, "n")
	if (is.null(shpfile)) {
		tmp <- attr(sn, "GeoDa")$shpfile
		if (is.null(tmp)) shpfile <- "unknown"
		else shpfile <- tmp
	}
	if (is.null(ind)) {
		tmp <- attr(sn, "GeoDa")$ind
		if (is.null(tmp)) ind <- "unknown"
		else ind <- tmp
	}
	con <- file(file, open="w")
	writeLines(paste("0", n, shpfile, ind, sep=" "), con)
	write.table(as.data.frame(sn), file=con, append=TRUE,
		row.names=FALSE, col.names=FALSE)
	close(con)
}
