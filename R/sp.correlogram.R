# Copyright 2002 by Roger Bivand
#

sp.correlogram <- function(neighbours, var, order=1, method="corr",
	style="W", zero.policy=FALSE) {
	if (class(neighbours) != "nb") stop("not a neighbours list")
	if (any(is.na(var))) stop("no NAs permitted in variable")
	nblags <- nblag(neighbours, maxlag=order)
	cardnos <- vector(mode="list", length=order)
	for (i in 1:order) cardnos[[i]] <- table(card(nblags[[i]]))
	if (method == "corr") {
		lags.x <- matrix(0, nrow=length(var), ncol=order)
		for (i in 1:order) lags.x[,i] <- lag.listw(nb2listw(nblags[[i]],
			style=style, zero.policy=zero.policy), var, 
			zero.policy=zero.policy)
		res <- cor(cbind(var, lags.x))[1,-1]
	} else if (method == "I") {
		res <- numeric(length=order)
		for (i in 1:order) {
			listw <- nb2listw(nblags[[i]], style=style, 
				zero.policy=zero.policy) 
			S0 <- Szero(listw)
			res[i] <- moran(var, listw, length(var), S0, 
				zero.policy=zero.policy)$I
		}			
	} else stop("method unknown")
	class(res) <- "spcor"
	names(res) <- 1:order
	attr(res, "method") <- method
	attr(res, "cardnos") <- cardnos
	attr(res, "var") <- deparse(substitute(var))
	res
}

print.spcor <- function(x, ...) {
	res <- as.vector(x)
	names(res) <- names(x)
	if (attr(x, "method") == "I") {
		meth <- "Moran's I"
	} else {
		meth <- "Spatial autocorrelation"
	}
	cat("Spatial correlogram for", attr(x, "var"), "\nmethod:", meth, "\n")
	print(res)
	invisible(res)
}

plot.spcor <- function(x, main, ylab, ylim, ...) {
	res <- as.vector(x)
	lags <- as.integer(names(x))
	if (missing(main)) main <- attr(x, "var")
	if (attr(x, "method") == "I") {
		if (missing(ylim)) ylim <- range(res)
		if (missing(ylab)) ylab <- "Moran's I"
	} else {
		if (missing(ylim)) ylim <- c(-1,1)
		if (missing(ylab)) ylab <- "Spatial autocorrelation"
	}
	plot(lags, res, type="h", ylim=ylim, main=main, ylab=ylab, 
		lwd=4, xaxt="n")
	axis(1, at=lags)
	abline(h=0)
}

