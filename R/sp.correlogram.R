# Copyright 2002-3 by Roger Bivand
#

sp.correlogram <- function (neighbours, var, order = 1, method = "corr", 
    style = "W", randomisation = TRUE, zero.policy = FALSE, spChk = NULL) {
    if (class(neighbours) != "nb") 
        stop("not a neighbours list")
    if (any(is.na(var))) 
        stop("no NAs permitted in variable")
    if (is.null(spChk)) 
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(var, nb2listw(neighbours, zero.policy = zero.policy))) 
        stop("Check of data and weights ID integrity failed")
    nblags <- nblag(neighbours, maxlag = order)
    cardnos <- vector(mode = "list", length = order)
    for (i in 1:order) cardnos[[i]] <- table(card(nblags[[i]]))
    if (method == "corr") {
        lags.x <- matrix(0, nrow = length(var), ncol = order)
        for (i in 1:order) lags.x[, i] <- lag.listw(nb2listw(nblags[[i]], 
            style = style, zero.policy = zero.policy), var, zero.policy = zero.policy)
        res <- cor(cbind(var, lags.x))[1, -1]
        names(res) <- 1:order
    }
    else if (method == "I") {
        res <- matrix(NA, nrow = order, ncol = 3)
        for (i in 1:order) {
            listw <- nb2listw(nblags[[i]], style = style, zero.policy = zero.policy)
            res[i,] <- moran.test(var, listw, randomisation = randomisation, 
		zero.policy = zero.policy)$estimate
        }
        rownames(res) <- 1:order
    }
    else stop("method unknown")
    obj <- list(res=res, method=method, cardnos=cardnos, var=deparse(substitute(var)))
    class(obj) <- "spcor"
    obj
}


#function(neighbours, var, order=1, method="corr",
#	style="W", zero.policy=FALSE, spChk=NULL) {
#	if (class(neighbours) != "nb") stop("not a neighbours list")
#	if (any(is.na(var))) stop("no NAs permitted in variable")
#	if (is.null(spChk)) spChk <- get.spChkOption()
#	if (spChk && !chkIDs(var, nb2listw(neighbours, 
#		zero.policy=zero.policy)))
#		stop("Check of data and weights ID integrity failed")
#	nblags <- nblag(neighbours, maxlag=order)
#	cardnos <- vector(mode="list", length=order)
#	for (i in 1:order) cardnos[[i]] <- table(card(nblags[[i]]))
#	if (method == "corr") {
#		lags.x <- matrix(0, nrow=length(var), ncol=order)
#		for (i in 1:order) lags.x[,i] <- lag.listw(nb2listw(nblags[[i]],
#			style=style, zero.policy=zero.policy), var, 
#			zero.policy=zero.policy)
#		res <- cor(cbind(var, lags.x))[1,-1]
#	} else if (method == "I") {
#		res <- numeric(length=order)
#		for (i in 1:order) {
#			listw <- nb2listw(nblags[[i]], style=style, 
#				zero.policy=zero.policy) 
#			S0 <- Szero(listw)
#			res[i] <- moran(var, listw, length(var), S0, 
#				zero.policy=zero.policy)$I
#		}			
#	} else stop("method unknown")
#	class(res) <- "spcor"
#	names(res) <- 1:order
#	attr(res, "method") <- method
#	attr(res, "cardnos") <- cardnos
#	attr(res, "var") <- deparse(substitute(var))
#	res
#}

print.spcor <- function (x, ...) 
{
    if (x$method == "I") {
        meth <- "Moran's I"
        res <- as.matrix(x$res)
        rownames(res) <- rownames(x$res)
        colnames(res) <- c("estimate", "expectation", "variance")
    } else {
        meth <- "Spatial autocorrelation"
        res <- as.vector(x$res)
        names(res) <- names(x$res)
    }
    cat("Spatial correlogram for", x$var, "\nmethod:", 
        meth, "\n")
    print(res)
    invisible(res)
}


#function(x, ...) {
#	res <- as.vector(x)
#	names(res) <- names(x)
#	if (attr(x, "method") == "I") {
#		meth <- "Moran's I"
#	} else {
#		meth <- "Spatial autocorrelation"
#	}
#	cat("Spatial correlogram for", attr(x, "var"), "\nmethod:", meth, "\n")
#	print(res)
#	invisible(res)
#}

plot.spcor <- function (x, main, ylab, ylim, ...) 
{
    if (missing(main)) 
        main <- x$var
    if (x$method == "I") {
        lags <- as.integer(rownames(x$res))
        if (missing(ylim)) 
            sd2 <- sqrt(x$res[,3])
            ylim <- range(c(x$res[,1]+sd2, x$res[,1]-sd2))
        if (missing(ylab)) 
            ylab <- "Moran's I"
        plot(lags, x$res[,1], type="p", pch=18, ylim = ylim, main = main, ylab = ylab, xaxt = "n")
#        segments(lags, x$res[,1], lags, x$res[,2], lwd=4, col="grey")
        arrows(lags, x$res[,1]+sd2, lags, x$res[,1]-sd2, length=0.1, angle=90)
        arrows(lags, x$res[,1]-sd2, lags, x$res[,1]+sd2, length=0.1, angle=90)
        axis(1, at = lags)
        abline(h = x$res[1,2])
    }
    else {
        res <- as.vector(x$res)
        lags <- as.integer(names(x$res))
        if (missing(ylim)) 
            ylim <- c(-1, 1)
        if (missing(ylab)) 
            ylab <- "Spatial autocorrelation"
        plot(lags, res, type = "h", ylim = ylim, main = main, ylab = ylab, 
            lwd = 4, xaxt = "n")
        axis(1, at = lags)
        abline(h = 0)
    }
}

#function(x, main, ylab, ylim, ...) {
#	res <- as.vector(x)
#	lags <- as.integer(names(x))
#	if (missing(main)) main <- attr(x, "var")
#	if (attr(x, "method") == "I") {
#		if (missing(ylim)) ylim <- range(res)
#		if (missing(ylab)) ylab <- "Moran's I"
#	} else {
#		if (missing(ylim)) ylim <- c(-1,1)
#		if (missing(ylab)) ylab <- "Spatial autocorrelation"
#	}
#	plot(lags, res, type="h", ylim=ylim, main=main, ylab=ylab, 
#		lwd=4, xaxt="n")
#	axis(1, at=lags)
#	abline(h=0)
#}

