# Copyright 1998-2004 by Roger Bivand and Andrew Bernat
#

lagsarlm <- function(formula, data = list(), listw, 
	na.action=na.fail, type="lag", method="eigen", quiet=TRUE, 
	zero.policy=FALSE, tol.solve=1.0e-7, tol.opt=.Machine$double.eps^0.5, 
	sparsedebug=FALSE) {
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, 
		method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	switch(type, lag = if (!quiet) cat("\nSpatial lag model\n"),
	    mixed = if (!quiet) cat("\nSpatial mixed autoregressive model\n"),
	    stop("\nUnknown model type\n"))
	if (!quiet) cat("Jacobian calculated using ")
	switch(method, 
		eigen = if (!quiet) cat("neighbourhood matrix eigenvalues\n"),
		sparse = if (!quiet) cat("sparse matrix techniques\n"),
		stop("...\nUnknown method\n"))
	y <- model.response(mf, "numeric")
#	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
#	if (any(is.na(x))) stop("NAs in independent variable")
	if (NROW(x) != length(listw$neighbours))
		stop("Input data and weights have different dimensions")
	n <- NROW(x)
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	if (any(is.na(wy))) stop("NAs in lagged dependent variable")
	if (type != "lag") {
		# check if there are enough regressors
	        if (m > 1) {
			WX <- matrix(nrow=n,ncol=(m-(K-1)))
			for (k in K:m) {
				wx <- lag.listw(listw, x[,k], 
				    zero.policy=zero.policy)
				if (any(is.na(wx))) 
				    stop("NAs in lagged independent variable")
				WX[,(k-(K-1))] <- wx
			}
		}
		if (K == 2) {
         	    # unnormalized weight matrices
                	if (!(listw$style == "W")) {
 	      			intercept <- as.double(rep(1, n))
       	        		wx <- lag.listw(listw, intercept, 
					zero.policy = zero.policy)
                    		if (m > 1) {
                        		WX <- cbind(wx, WX)
                    		} else {
			      		WX <- matrix(wx, nrow = n, ncol = 1)
                    		}
                	} 
            	}   
		m1 <- m + 1
		mm <- NCOL(x) + NCOL(WX)
            	xxcolnames <- character(mm)
		for (k in 1:m) xxcolnames[k] <- xcolnames[k]
		for (k in m1:mm) 
		    xxcolnames[k] <- paste("lag.", xcolnames[k-mm+m], sep="")
		x <- cbind(x, WX)
		colnames(x) <- xxcolnames
		m <- NCOL(x)
		rm(wx, WX)
	}
	if (method == "eigen") {
		if (!quiet) cat("Computing eigenvalues ...\n")
		eig <- eigenw(listw)
		if (!quiet) cat("\n")
#range inverted 031031, email from Salvati Nicola (and Rein Halbersma)
		if (is.complex(eig)) eig.range <- 1/range(Re(eig))
		else eig.range <- 1/range(eig)
		lm.null <- lm(y ~ x - 1)
		lm.w <- lm.fit(x, wy)
		e.null <- lm.null$residuals
		e.w <- lm.w$residuals
		e.a <- t(e.null) %*% e.null
		e.b <- t(e.w) %*% e.null
		e.c <- t(e.w) %*% e.w
		opt <- optimize(sar.lag.mixed.f, interval=eig.range,
			maximum=TRUE, tol=tol.opt, eig=eig,
			e.a=e.a, e.b=e.b, e.c=e.c, n=n, quiet=quiet)
	} else {
		opt <- dosparse(listw, y, x, wy, K, quiet, tol.opt, 
			sparsedebug)
	}
	rho <- c(opt$maximum)
	names(rho) <- "rho"
	LL <- c(opt$objective)
	lm.lag <- lm((y - rho*wy) ~ x - 1)
	r <- residuals(lm.lag)
	fit <- y - r
	names(r) <- names(fit)
	coef.rho <- coefficients(lm.lag)
	names(coef.rho) <- colnames(x)
	SSE <- deviance(lm.lag)
	s2 <- SSE/n
	if (method != "eigen") {
		LLs <- opt$LLs
		lm.null <- opt$lm.null
		rest.se <- NULL
		rho.se <- NULL
		LMtest <- NULL
		ase <- FALSE
		varb <- FALSE
	} else {
		LLs <- NULL
		tr <- function(A) sum(diag(A))
# beware of complex eigenvalues!
		O <- (eig/(1-rho*eig))^2
		omega <- sum(O)
		if (is.complex(omega)) omega <- Re(omega)
		W <- listw2mat(listw)
		A <- solve(diag(n) - rho*W, tol=tol.solve)
		AW <- A %*% W
		zero <- rbind(rep(0,length(coef.rho)))
		xtawxb <- s2*(t(x) %*% AW %*% x %*% coef.rho)
		V <- s2*(s2*tr(t(AW) %*% AW) +
			t(AW %*% x %*% coef.rho) %*%
			(AW %*% x %*% coef.rho)) + omega*s2^2
		inf1 <- rbind(n/2, s2*tr(AW), t(zero))
		inf2 <- rbind(s2*tr(AW), V, xtawxb)
		xtx <- s2*t(x) %*% x
		inf3 <- rbind(zero, t(xtawxb), xtx)
		inf <- cbind(inf1, inf2, inf3)
		varb <- (s2^2) * solve(inf, tol=tol.solve)
		rownames(varb) <- colnames(varb) <- 
			c("sigma", "rho", colnames(x))
		rest.se <- sqrt(diag(varb))[-c(1:2)]
		rho.se <- sqrt(varb[2,2])
		TW <- (W %*% W) + (t(W) %*% W)
		T22 <- sum(diag(TW))
		T21A <- sum(diag(TW %*% A))
		LMtest <- ((t(r) %*% W %*% r)/s2)^2
		LMtest <- LMtest/(T22 - ((T21A^2)*(rho.se^2)))
		ase <- TRUE
	}
	call <- match.call()
	ret <- structure(list(type=type, rho=rho, 
		coefficients=coef.rho, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), lm.model=lm.null,
		method=method, call=call, residuals=r, 
		lm.target=lm.lag, fitted.values=fit,
		se.fit=NULL, formula=formula,
		ase=ase, LLs=LLs, rho.se=rho.se, LMtest=LMtest, 
		resvar=varb, zero.policy=zero.policy), class=c("sarlm"))
	if (zero.policy) {
		zero.regs <- attr(listw$neighbours, 
			"region.id")[which(card(listw$neighbours) == 0)]
		if (length(zero.regs) > 0)
			attr(ret, "zero.regs") <- zero.regs
	}
	if (!is.null(na.act))
		ret$na.action <- na.act
	ret
}

sar.lag.mixed.f <- function(rho, eig, e.a, e.b, e.c, n, quiet)
{
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
	s2 <- SSE/n
	if (is.complex(eig)) det <- Re(prod(1 - rho*eig)) 
	else det <- prod(1 - rho*eig)
	ret <- (log(det) - ((n/2)*log(2*pi)) - (n/2)*log(s2)
		- (1/(2*s2))*SSE)
	if (!quiet) cat("Rho:\t", rho, "\tfunction value:\t", ret, "\n")
	ret
}

sar.lag.mixed.f.s <- function(rho, sn, e.a, e.b, e.c, n, quiet, sparsedebug)
{
	SSE <- e.a - 2*rho*e.b + rho*rho*e.c
	s2 <- SSE/n
	ret <- (logSpwdet(sparseweights=sn, rho=rho, debug=sparsedebug)
		- ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
	if (!quiet) cat("Rho:\t", rho, "\tfunction value:\t", ret, "\n")
	ret
}


dosparse <- function (listw, y, x, wy, K, quiet, tol.opt, sparsedebug) {
	sn <- listw2sn(listw)
	m <- ncol(x)
	n <- nrow(x)
	LLs <- vector(mode="list", length=length(K:m))
	j <- 1
	for (i in K:m) {
		thisx <- x[,-i]
		lm.null <- lm.fit(thisx, y)
		lm.w <- lm.fit(thisx, wy)
		e.null <- lm.null$residuals
		e.w <- lm.w$residuals
		e.a <- t(e.null) %*% e.null
		e.b <- t(e.w) %*% e.null
		e.c <- t(e.w) %*% e.w
		LLs[[j]] <- optimize(sar.lag.mixed.f.s, interval=c(-1,1),
		maximum=TRUE, tol=tol.opt, sn=sn,
		e.a=e.a, e.b=e.b, e.c=e.c, n=n, quiet=quiet,
		sparsedebug=sparsedebug)$objective
		attr(LLs[[j]], "nall") <- n
		attr(LLs[[j]], "nobs") <- n
		attr(LLs[[j]], "df") <- (m+2)-1
		attr(LLs[[j]], "name") <- colnames(x)[i]
		class(LLs[[j]]) <- "logLik"
		j <- j + 1
	}
	lm.null <- lm(y ~ x - 1)
	lm.w <- lm.fit(x, wy)
	e.null <- lm.null$residuals
	e.w <- lm.w$residuals
	e.a <- t(e.null) %*% e.null
	e.b <- t(e.w) %*% e.null
	e.c <- t(e.w) %*% e.w
	sn <- listw2sn(listw)
	opt <- optimize(sar.lag.mixed.f.s, interval=c(-1,1),
		maximum=TRUE, tol=tol.opt, sn=sn,
		e.a=e.a, e.b=e.b, e.c=e.c, n=n, quiet=quiet, 
		sparsedebug=sparsedebug)
	maximum <- opt$maximum
	objective <- opt$objective
	res <- list(maximum=maximum, objective=objective, LLs=LLs,
		lm.null=lm.null)
}
