# Copyright 2001-2003 by Roger Bivand with contributions by Stéphane Dray
#
	


poly2nb <- function(pl, row.names=NULL, snap=sqrt(.Machine$double.eps),
	queen=TRUE) {
	if (!inherits(pl, "polylist")) stop("Not a polygon list")
	if (inherits(pl, "multiparts")) stop("Convert to newer polylist format")
	n <- length(pl)
	regid <- attr(pl, "region.id")
	if (is.null(regid)) {
		if(is.null(row.names)) regid <- as.character(1:n)
		else {
			if(length(row.names) != n)
				stop("row.names wrong length")
			else if (length(unique(row.names)) != length(row.names))
	    			stop("non-unique row.names given")
			else regid <- row.names
		}
	}
	poly2bbs <- function(pl) {
		n <- length(pl)
		res <- matrix(0, nrow=n, ncol=4)
		for (i in 1:n) res[i,] <- attr(pl[[i]], "bbox")
		res
	}
	bb <- poly2bbs(pl)
	

	between <- function(x, low, up) {return(x >= low && x <= up)}

	pipbb <- function(pt, bb) {
		return(between(pt[1], bb[1], bb[3]) && 
			between(pt[2], bb[2], bb[4]))
	}

	polypoly <- function(poly1, poly2, snap) {
		snap2 <- snap^2
		a <- outer(poly1[,1], poly2[,1], "-")
		b <- outer(poly1[,2], poly2[,2], "-")
		c <- a^2 + b^2
		d <- rle(unlist(apply(c,2,function(x) {
			length(which(x < snap2)) > 0} )))
		res <- length(which(d$values == TRUE))
		if (res > 0) res <- sum(d$length[d$values == TRUE])
		res
	}

	pipbbij <- function(bbi, bbj) {
		jhit <- logical(8)
	    	jhit[1] <- pipbb(c(bbi[1], bbi[2]), bbj)
		jhit[2] <- pipbb(c(bbi[1], bbi[4]), bbj)
		jhit[3] <- pipbb(c(bbi[3], bbi[2]), bbj)
		jhit[4] <- pipbb(c(bbi[3], bbi[4]), bbj)
		jhit[5] <- pipbb(c(bbj[1], bbj[2]), bbi)
		jhit[6] <- pipbb(c(bbj[1], bbj[4]), bbi)
		jhit[7] <- pipbb(c(bbj[3], bbj[2]), bbi)
		jhit[8] <- pipbb(c(bbj[3], bbj[4]), bbi)
		return(jhit)
	}

	ans <- vector(mode="list", length=n)
	for (i in 1:n) ans[[i]] <- integer(0)
	criterion <- ifelse(queen, 0, 1)
	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			jhit <- pipbbij(bb[i,], bb[j,])
			if (any(jhit)) {
			    khit <- 0
			    khit <- polypoly(na.omit(pl[[i]]), 
				na.omit(pl[[j]]), snap)
			    if (khit > criterion) {
				ans[[i]] <- c(ans[[i]], j)
				ans[[j]] <- c(ans[[j]], i)
			    }
			}
		}
	}
	for (i in 1:n) ans[[i]] <- sort(ans[[i]])
	class(ans) <- "nb"
	attr(ans, "region.id") <- regid
	attr(ans, "call") <- match.call()
	if (queen) attr(ans, "type") <- "queen"
	else attr(ans, "type") <- "rook"
	invisible(ans)
}

#shape2poly <- function(shape, region.id=NULL) {
#    if (is.null(shape$shp)) stop("No shp component in this list")
#    if (shape$shp$header$shape.type != 5) stop("Not a polygon shapefile")
#    nrecord <- length(shape$shp$shp)
#    res <- vector(mode="list", length=nrecord)
#    if (!is.null(region.id)) {
#	if (length(unique(region.id)) == nrecord) id <- region.id
#	else region.id <- NULL
#    }
#    if (is.null(region.id)) id <- vector(mode="character", length=nrecord)
#    for (i in 1:nrecord) {
#	res[[i]] <- as.matrix(shape$shp$shp[[i]]$points)
#	if (is.null(region.id)) id[i]<- as.character(shape$dbf$dbf[i,1])
#    }
#
#    attr(res, "region.id") <- id
#    class(res) <- "polylist"
#    return(res)
#
#}
#
#shape2bbs <- function(shape) {
#    if (is.null(shape$shp)) stop("No shp component in this list")
#    if (shape$shp$header$shape.type != 5) stop("Not a polygon shapefile")
#    n <- length(shape$shp$shp)
#    res <- matrix(0, ncol=4, nrow=n)
#    for (i in 1:n) res[i,] <- as.vector(shape$shp$shp[[i]]$box)
#    res
#}
#
#Map2poly <- function(Map, region.id=NULL) {
#	res <- .get.polylist(Map=Map, region.id=region.id)
#	res
#}
#
#.get.polylist <- function(Map, region.id=NULL) {
#	if (class(Map) != "Map") stop("not a Map")
#	n <- length(Map$Shapes)
#	res <- vector(mode="list", length=n)
#	nParts <- integer(n)
#	for (i in 1:n) nParts[i] <- attr(Map$Shapes[[i]], "nParts")
#	if (any(nParts != 1)) {
#		for (i in 1:n) {
#			Pstart <- Map$Shapes[[i]]$Pstart
#			nVerts <- attr(Map$Shapes[[i]], "nVerts")
#			from <- integer(nParts[i])
#			to <- integer(nParts[i])
#			from[1] <- 1
#			for (j in 1:nParts[i]) {
#				if (j == nParts[i]) to[j] <- nVerts
#				else {
#					to[j] <- Pstart[j+1]
#					from[j+1] <- to[j]+1
#				}
#			}
#			res[[i]] <- Map$Shapes[[i]]$verts[from[1]:to[1],]
#			for (j in 2:nParts[i]) {
#			    res[[i]] <- rbind(res[[i]], c(NA, NA))
#			    res[[i]] <- rbind(res[[i]], 
#				Map$Shapes[[i]]$verts[from[j]:to[j],])
#			}
#		}
#	} else {
#		for (i in 1:n) res[[i]] <- Map$Shapes[[i]]$verts
#	}
#	if (is.null(region.id) || length(region.id) != n) {
#		attr(res, "region.id") <- as.character(1:n)
#	} else {
#		attr(res, "region.id") <- as.character(region.id)
#	}
#	class(res) <- "polylist"
#	invisible(res)
#}
#
#convert.pl <- function(pl) {
#	if (!inherits(pl, "multiparts")) stop("not a mulitpart polylist")
#	res <- vector(mode="list", length=length(pl))
#	for (i in 1:length(pl)) {
#		lp <- length(pl[[i]])
#		res[[i]] <- pl[[i]][[1]]
#		if (lp > 1) {
#			for (j in 2:lp) {
#				res[[i]] <- rbind(res[[i]], c(NA, NA))
#				res[[i]] <- rbind(res[[i]], pl[[i]][[j]])
#			}
#		}
#	}
#	if (!is.null(attr(pl, "region.id")))
#		attr(res, "region.id") <- attr(pl, "region.id")
#	class(res) <- "polylist"
#	res
#}
#
#.get.polybbs <- function(Map) {
#	if (class(Map) != "Map") stop("not a Map")
#	n <- length(Map$Shapes)
#	res <- matrix(0, ncol=4, nrow=n)
#	for (i in 1:n) res[i,] <- attr(Map$Shapes[[i]], "bbox")
#	res
#}
#
#Map2bbs <- function(Map) {
#	res <- .get.polybbs(Map)
#	res
#}
#
#leglabs <- function(vec, under="under", over="over", between="-") {
#	x <- vec
#	res <- character(length(x)-1)
#	res[1] <- paste(under, x[2])
#	for (i in 2:(length(x)-2)) res[i] <- paste(x[i], between, x[i+1])
#	res[length(x)-1] <- paste(over, x[length(x)-1])
#	res
#}
#
#findInterval2 <- function (x, vec, rightmost.closed = FALSE, all.inside = TRUE) 
#{
#    nx <- length(x)
#    if (any(is.na(vec) | is.nan(vec))) stop ("NAs found in vec")
#    if (is.unsorted(vec)) 
#        stop("`vec' must be sorted non-decreasingly")
#    if (vec[1] == -Inf) vec[1] <- -(.Machine$double.xmax)
#    if (vec[length(vec)] == Inf) 
#	vec[length(vec)] <- .Machine$double.xmax
#    .C("find_interv_vec", xt = as.double(vec), n = length(vec), 
#        x = as.double(x), nx = nx, as.logical(rightmost.closed), 
#        as.logical(all.inside), index = integer(nx), DUP = FALSE,
#	PACKAGE = "base")$index
#}
#
