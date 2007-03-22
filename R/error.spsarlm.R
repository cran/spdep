# Copyright 1998-2007 by Roger Bivand (non-W styles Rein Halbersma)
#

errorsarlm <- function(formula, data = list(), listw, na.action=na.fail, 
	method="eigen", quiet=TRUE, zero.policy=FALSE, interval=c(-1,0.999), 
	tol.solve=1.0e-10, tol.opt=.Machine$double.eps^0.5, cholAlloc=NULL) {
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	can.sim <- as.logical(NA)
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	if (!quiet) cat(paste("\nSpatial autoregressive error model\n", 
		"Jacobian calculated using "))
	switch(method,
		eigen = if (!quiet) cat("neighbourhood matrix eigenvalues\n"),
	        SparseM = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		    stop("SparseM method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		    stop("SparseM method requires symmetric weights")
		    if (!quiet) cat("sparse matrix techniques using SparseM\n")
		},
	        Matrix = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		    stop("Matrix method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		    stop("Matrix method requires symmetric weights")
		    if (!quiet) cat("sparse matrix techniques using Matrix\n")
		},
		stop("...\n\nUnknown method\n"))
	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	n <- NROW(x)
	m <- NCOL(x)
# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
	if (m > 1) {
	    WX <- matrix(nrow=n,ncol=(m-(K-1)))
	    for (k in K:m) {
		wx <- lag.listw(listw, x[,k], zero.policy=zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
	    }
	}
	if (K == 2) {
# modified to meet other styles, email from Rein Halbersma
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
		if (m > 1) WX <- cbind(wx, WX)
		else WX <- matrix(wx, nrow=n, ncol=1)
	}
	colnames(WX) <- xcolnames
	rm(wx)
	similar <- FALSE
	if (method == "eigen") {
		if (!quiet) cat("Computing eigenvalues ...\n")
		if (listw$style %in% c("W", "S") & can.sim) {
			eig <- eigenw(similar.listw(listw))
			similar <- TRUE
		} else eig <- eigenw(listw)
		if (!quiet) cat("\n")
# range inverted 031031, email from Salvati Nicola (and Rein Halbersma)
		if (is.complex(eig)) eig.range <- 1/range(Re(eig))
		else eig.range <- 1/range(eig)
		opt <- optimize(sar.error.f, 
			lower=eig.range[1]+.Machine$double.eps, 
			upper=eig.range[2]-.Machine$double.eps, maximum=TRUE,
			tol=tol.opt, eig=eig, y=y, wy=wy, x=x, WX=WX, 
			n=n, quiet=quiet)
		lambda <- opt$maximum
		names(lambda) <- "lambda"
		LL <- opt$objective
	} else if (method == "SparseM") {
		if (listw$style %in% c("W", "S") & can.sim) {
#			csrw <- asMatrixCsrListw(similar.listw(listw))
			csrw <- asMatrixCsrListw(similar.listw(listw),
        			zero.policy=zero.policy)
			similar <- TRUE
		} else csrw <- asMatrixCsrListw(listw, zero.policy=zero.policy)
#		} else csrw <- asMatrixCsrListw(listw)
		gc(FALSE)
		I <- asMatrixCsrI(n)
		if (is.null(cholAlloc)) {
		# Martin Reismann large sparse nnzlmax problem
			nlink <- sum(card(listw$neighbours))
			tmpmax <- 3 * (nlink + n)
			nnzlmax <- max(10*nlink, floor(.2*nlink^1.4))
			nsubmax <- tmpmax
			cholAlloc <- list(nsubmax=nsubmax, nnzlmax=nnzlmax,
				tmpmax=tmpmax)
		}
		# tmpmax and gc() calls: Danlin Yu 20041213
		opt <- optimize(sar.error.f.sM, interval=interval, 
			maximum=TRUE, tol=tol.opt, csrw=csrw, I=I, y=y, wy=wy, 
			x=x, WX=WX, n=n, cholAlloc=cholAlloc, quiet=quiet)
		lambda <- opt$maximum
		names(lambda) <- "lambda"
		LL <- opt$objective
		gc(FALSE)
	} else if (method == "Matrix") {
        	if (listw$style %in% c("W", "S") & can.sim) {
	    	    csrw <- as_dsTMatrix_listw(listw2U(similar.listw(listw)))
	    	    similar <- TRUE
		} else csrw <- as_dsTMatrix_listw(listw)
		gc(FALSE)
        	I <- as_dgCMatrix_I(n)
		I <- as(I, "CsparseMatrix")
		opt <- optimize(sar.error.f.M, interval=interval, 
			maximum=TRUE, tol=tol.opt, csrw=csrw, I=I, y=y, wy=wy, 
			x=x, WX=WX, n=n, quiet=quiet)
		lambda <- opt$maximum
		names(lambda) <- "lambda"
		LL <- opt$objective
		gc(FALSE)
	}
	lm.target <- lm(I(y - lambda*wy) ~ I(x - lambda*WX) - 1)
	r <- as.vector(residuals(lm.target))
	fit <- as.vector(y - r)
	p <- lm.target$rank
	SSE <- deviance(lm.target)
	s2 <- SSE/n
	rest.se <- (summary(lm.target)$coefficients[,2])*sqrt((n-p)/n)
	coef.lambda <- coefficients(lm.target)
	names(coef.lambda) <- xcolnames
	lm.model <- lm(formula, data)
	ase <- FALSE
	lambda.se <- NULL
	LMtest <- NULL
	asyvar1 <- FALSE
	if (method == "eigen") {
		tr <- function(A) sum(diag(A))
		W <- listw2mat(listw)
		A <- solve(diag(n) - lambda*W)
		WA <- W %*% A
		asyvar <- matrix(0, nrow=2+p, ncol=2+p)
		asyvar[1,1] <- n / (2*(s2^2))
		asyvar[2,1] <- asyvar[1,2] <- tr(WA) / s2
		asyvar[2,2] <- tr(WA %*% WA) + tr(t(WA) %*% WA)
		asyvar[3:(p+2),3:(p+2)] <- s2*(t(x - lambda*WX) %*% 
			(x - lambda*WX))
		asyvar1 <- solve(asyvar, tol=tol.solve)
		rownames(asyvar1) <- colnames(asyvar1) <- 
			c("sigma", "lambda", xcolnames)
		
		lambda.se <- sqrt(asyvar1[2,2])
		ase <- TRUE
	}
	call <- match.call()
	names(r) <- names(y)
	names(fit) <- names(y)
	ret <- structure(list(type="error", lambda=lambda,
		coefficients=coef.lambda, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), lm.model=lm.model, 
		method=method, call=call, residuals=r, lm.target=lm.target,
		opt=opt, fitted.values=fit, ase=ase, formula=formula,
		se.fit=NULL, resvar=asyvar1, similar=similar,
		lambda.se=lambda.se, LMtest=LMtest, zero.policy=zero.policy, 
		aliased=aliased), class=c("sarlm"))
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

sar.error.f <- function(lambda, eig, y, wy, x, WX, n, quiet)
{
	yl <- y - lambda*wy
	xl <- x - lambda*WX
	xl.q <- qr.Q(qr(xl))
	xl.q.yl <- t(xl.q) %*% yl
	SSE <- t(yl) %*% yl - t(xl.q.yl) %*% xl.q.yl
	s2 <- SSE/n
	if (is.complex(eig)) det <- Re(prod(1 - lambda*eig)) 
	else det <- prod(1 - lambda*eig)
	ret <- (log(det) - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
	if (!quiet) cat("lambda:", lambda, " function:", ret, " Jacobian:", log(det), " SSE:", SSE, "\n")
	ret
}


sar.error.f.sM <- function(lambda, csrw, I, y, wy, x, WX, n, cholAlloc, quiet) {
	yl <- y - lambda*wy
	xl <- x - lambda*WX
	xl.q <- qr.Q(qr(xl))
	xl.q.yl <- t(xl.q) %*% yl
	SSE <- t(yl) %*% yl - t(xl.q.yl) %*% xl.q.yl
	s2 <- SSE/n
        Det <- get("det", "package:SparseM")
	Jacobian <- log(Det(chol((I - lambda * csrw), 
		nsubmax=cholAlloc$nsubmax, nnzlmax=cholAlloc$nnzlmax, 
		tmpmax=cholAlloc$tmpmax))^2)
	gc(FALSE)
	ret <- (Jacobian -
		((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
	if (!quiet) cat("lambda:", lambda, " function:", ret, " Jacobian:", Jacobian, " SSE:", SSE, "\n")
	ret
}

sar.error.f.M <- function(lambda, csrw, I, y, wy, x, WX, n, quiet) {
	yl <- y - lambda*wy
	xl <- x - lambda*WX
	xl.q <- qr.Q(qr(xl))
	xl.q.yl <- t(xl.q) %*% yl
	SSE <- t(yl) %*% yl - t(xl.q.yl) %*% xl.q.yl
	s2 <- SSE/n
    	CHOL <- try(chol(as((I - lambda * csrw), "dsCMatrix")), silent=TRUE)
    	if (class(CHOL) == "try-error") {
        	Jacobian <- NA
    	} else {
        	Jacobian <- sum(2*log(diag(CHOL)))
    	}
	gc(FALSE)
	ret <- (Jacobian -
		((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
	if (!quiet) cat("lambda:", lambda, " function:", ret, " Jacobian:", Jacobian, " SSE:", SSE, "\n")
	ret
}

