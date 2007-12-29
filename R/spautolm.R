# Copyright 2005-7 by Roger Bivand
spautolm <- function(formula, data = list(), listw, weights,
    na.action=na.fail, verbose=FALSE, tol.opt=.Machine$double.eps^(2/3),
    family="SAR", method="full", interval=c(-1,0.999), zero.policy=FALSE,
#    cholAlloc=NULL, 
    tol.solve=.Machine$double.eps, llprof=NULL) 
{
    if (!inherits(listw, "listw")) 
        stop("No neighbourhood list")

    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")

#    mt <- terms(formula, data = data)
#    mf <- lm(formula, data, , weights, na.action=na.action,
#        method="model.frame")
    na.act <- attr(mf, "na.action")
    if (!is.null(na.act)) {
        subset <- !(1:length(listw$neighbours) %in% na.act)
        listw <- subset(listw, subset, zero.policy=zero.policy)
    }

    Y <- model.extract(mf, "response")
    if (any(is.na(Y))) stop("NAs in dependent variable")
    X <- model.matrix(mt, mf)
    if (any(is.na(X))) stop("NAs in independent variable")
    n <- nrow(X)
    weights <- as.vector(model.extract(mf, "weights"))
# set up default weights
    if (!is.null(weights) && !is.numeric(weights)) 
        stop("'weights' must be a numeric vector")
    if (is.null(weights)) weights <- rep(as.numeric(1), n)
    if (any(is.na(weights))) stop("NAs in weights")
    if (any(weights < 0)) stop("negative weights")
    lm.base <- lm(Y ~ X - 1, weights=weights)
    aliased <- is.na(coefficients(lm.base))
    cn <- names(aliased)
    names(aliased) <- substr(cn, 2, nchar(cn))
    if (any(aliased)) {
        nacoef <- which(aliased)
	x <- x[,-nacoef]
    }
    can.sim <- as.logical(NA)
    if (listw$style %in% c("W", "S")) 
	can.sim <- can.be.simmed(listw)

    sum_lw <- sum(log(weights))
    if (!is.null(llprof) && method != "full" && length(llprof) == 1)
        stop("sequence of values required for llprof")
    if (method == "full") {
# spatial weights matrix
        W <- listw2mat(listw)
	attr(W, "dimnames") <- NULL
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
# range for line search
        eig <- eigen(diag(sqrt(weights)) %*% W %*% diag(sqrt(1/weights)), 
            only.values =TRUE)$values
        if (is.complex(eig)) eig.range <- 1/range(Re(eig))
        else eig.range <- 1/range(eig)
        I <- diag(n)
# do line search
        dweights <- diag(weights)
        if (!is.null(llprof)) {
            if (length(llprof) == 1)
                llprof <- seq(eig.range[1]+.Machine$double.eps,
                    eig.range[2]-.Machine$double.eps, length.out=llprof)
            ll_prof <- numeric(length(llprof))
            for (i in seq(along=llprof)) ll_prof[i] <- .opt.fit.full(
                llprof[i], Y=Y, X=X, n=n, W=W, eig=eig, I=I,
                weights=dweights, sum_lw=sum_lw, family=family,
                verbose=verbose, tol.solve=tol.solve)
        }
        opt <- optimize(.opt.fit.full, lower=eig.range[1]+.Machine$double.eps,
            upper=eig.range[2]-.Machine$double.eps, maximum=TRUE,
            tol = tol.opt, Y=Y, X=X, n=n, W=W, eig=eig, I=I,
            weights=dweights, sum_lw=sum_lw, family=family,
            verbose=verbose, tol.solve=tol.solve)
        lambda <- opt$maximum
        names(lambda) <- "lambda"
        LL <- opt$objective
# get GLS coefficients
        fit <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, I=I,
            weights=diag(weights), family=family, out=TRUE,
		tol.solve=tol.solve)
# create residuals and fitted values (Cressie 1993, p. 564)
	fit$signal_trend <- drop(X %*% fit$coefficients)
	fit$signal_stochastic <- drop(lambda * W %*% (Y - fit$signal_trend))
	fit$fitted.values <- fit$signal_trend + fit$signal_stochastic
	fit$residuals <- drop(Y - fit$fitted.values)
# get null LL
        LL0 <- .opt.fit.full(lambda=0, Y=Y, X=X, n=n, W=W, eig=eig, I=I,
            weights=diag(weights), sum_lw=sum_lw, family=family, verbose=FALSE,
		tol.solve=tol.solve)
#    \} else if (method == "SparseM") \{
#        if (family == "SMA") stop("SMA only for full method")
#        if (listw$style %in% c("W", "S") && !can.sim)
#        stop("SparseM method requires symmetric weights")
#        if (listw$style %in% c("B", "C", "U") && 
# 	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
#	    stop("SparseM method requires symmetric weights")
#        I <- asMatrixCsrI(n)
#	W <- asMatrixCsrListw(listw, zero.policy=zero.policy)
#        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
#	    warning("Non-symmetric spatial weights in CAR model")
## Jacobian only from symmetric W_J, W can be asymmetric though
#        if (listw$style %in% c("W", "S") & can.sim) {
#	    W_J <- asMatrixCsrListw(similar.listw(listw),
#                zero.policy=zero.policy)
##	    similar <- TRUE
#	} else W_J <- W
#	gc(FALSE)
#        Sweights <- new("matrix.csr", ra=weights, ja=1:n, ia=1:(n+1), 
#	    dimension=c(n,n))
##	tmpmax <- sum(card(listw$neighbours)) + n
#	if (is.null(cholAlloc)) {
#	# Martin Reismann large sparse nnzlmax problem
#		nlink <- sum(card(listw$neighbours))
#		tmpmax <- 3 * (nlink + n)
#		nnzlmax <- max(10*nlink, floor(.2*nlink^1.4))
#		nsubmax <- tmpmax
#		cholAlloc <- list(nsubmax=nsubmax, nnzlmax=nnzlmax,
#			tmpmax=tmpmax)
#	}
## do line search
#        if (!is.null(llprof)) {
#            ll_prof <- numeric(length(llprof))
#            for (i in seq(along=llprof)) ll_prof[i] <- .opt.fit.SparseM(
#                llprof[i], Y=Y, X=X, n=n, W=W, W_J=W_J, I=I,
#                weights=Sweights, sum_lw=sum_lw, family=family,
#                verbose=verbose, cholAlloc=cholAlloc, tol.solve=tol.solve)
#        }
#        opt <- optimize(.opt.fit.SparseM, lower=interval[1],
#            upper=interval[2], maximum=TRUE,
#            tol = tol.opt, Y=Y, X=X, n=n, W=W, W_J=W_J, I=I,
#            weights=Sweights, sum_lw=sum_lw, family=family,
#            verbose=verbose, cholAlloc=cholAlloc, tol.solve=tol.solve)
#        lambda <- opt$maximum
#        names(lambda) <- "lambda"
#        LL <- opt$objective
## get GLS coefficients
#        fit <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, I=I,
#            weights=Sweights, family=family, out=TRUE, tol.solve=tol.solve)
## create residuals and fitted values (Cressie 1993, p. 564)
#	fit$signal_trend <- drop(X %*% fit$coefficients)
#	fit$signal_stochastic <- drop(lambda * W %*% (Y - fit$signal_trend))
#	fit$fitted.values <- fit$signal_trend + fit$signal_stochastic
#	fit$residuals <- drop(Y - fit$fitted.values)
## get null LL
#        LL0 <- .opt.fit.SparseM(lambda=as.numeric(0), Y=Y, X=X, n=n, W=W,
#            W_J=W_J, I=I, weights=Sweights, sum_lw=sum_lw, family=family,
#            verbose=verbose, cholAlloc=cholAlloc, tol.solve=tol.solve)
##        weights <- diag(Sweights)
    }  else if (method == "Matrix") {
        if (family == "SMA") stop("SMA only for full method")
        if (listw$style %in% c("W", "S") && !can.sim)
        stop("Matrix method requires symmetric weights")
        if (listw$style %in% c("B", "C", "U") && 
 	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
	    stop("Matrix method requires symmetric weights")
        I <- as_dgCMatrix_I(n)
	I <- as(I, "CsparseMatrix")
	W <- as_dgRMatrix_listw(listw)
	W <- as(W, "CsparseMatrix")
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
# Jacobian only from symmetric W_J, W can be asymmetric though
        if (listw$style %in% c("W", "S") & can.sim) {
	    W_J <- as_dsTMatrix_listw(listw2U(similar.listw(listw)))
#	    similar <- TRUE
	} else W_J <- as_dsTMatrix_listw(listw)
#	gc(FALSE)
        Sweights <- as(Diagonal(x=weights), "sparseMatrix")
# do line search
        if (!is.null(llprof)) {
            ll_prof <- numeric(length(llprof))
            for (i in seq(along=llprof)) ll_prof[i] <- .opt.fit.Matrix(
                llprof[i], Y=Y, X=X, n=n, W=W, W_J=W_J, I=I,
                weights=Sweights, sum_lw=sum_lw, family=family,
                verbose=verbose, tol.solve=tol.solve)
        }
        opt <- optimize(.opt.fit.Matrix, lower=interval[1],
            upper=interval[2], maximum=TRUE,
            tol = tol.opt, Y=Y, X=X, n=n, W=W, W_J=W_J, I=I,
            weights=Sweights, sum_lw=sum_lw, family=family,
            verbose=verbose, tol.solve=tol.solve)
        lambda <- opt$maximum
        names(lambda) <- "lambda"
        LL <- opt$objective
# get GLS coefficients
        fit <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, I=I,
            weights=Sweights, family=family, out=TRUE, tol.solve=tol.solve)
# create residuals and fitted values (Cressie 1993, p. 564)
	fit$signal_trend <- drop(X %*% fit$coefficients)
	fit$signal_stochastic <- drop(lambda * W %*% (Y - fit$signal_trend))
	fit$fitted.values <- fit$signal_trend + fit$signal_stochastic
	fit$residuals <- drop(Y - fit$fitted.values)
# get null LL
        LL0 <- .opt.fit.Matrix(lambda=as.numeric(0), Y=Y, X=X, n=n, W=W,
            W_J=W_J, I=I, weights=Sweights, sum_lw=sum_lw, family=family,
            verbose=verbose, tol.solve=tol.solve)
#        weights <- diag(Sweights)
    }  else if (method == "spam") {
        if (family == "SMA") stop("SMA only for full method")
        if (listw$style %in% c("W", "S") && !can.sim)
        stop("spam method requires symmetric weights")
        if (listw$style %in% c("B", "C", "U") && 
 	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
	    stop("spam method requires symmetric weights")
        I <- diag.spam(1, n, n)
	W <- as.spam.listw(listw)
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
# Jacobian only from symmetric W_J, W can be asymmetric though
        if (listw$style %in% c("W", "S") & can.sim) {
	    W_J <- as.spam.listw(listw2U(similar.listw(listw)))
#	    similar <- TRUE
	} else W_J <- W
#	gc(FALSE)
        Sweights <- diag.spam(x=weights, n, n)
# do line search
        if (!is.null(llprof)) {
            ll_prof <- numeric(length(llprof))
            for (i in seq(along=llprof)) ll_prof[i] <- .opt.fit.spam(
                llprof[i], Y=Y, X=X, n=n, W=W, W_J=W_J, I=I,
                weights=Sweights, sum_lw=sum_lw, family=family,
                verbose=verbose, tol.solve=tol.solve)
        }
        opt <- optimize(.opt.fit.spam, lower=interval[1],
            upper=interval[2], maximum=TRUE,
            tol = tol.opt, Y=Y, X=X, n=n, W=W, W_J=W_J, I=I,
            weights=Sweights, sum_lw=sum_lw, family=family,
            verbose=verbose, tol.solve=tol.solve)
        lambda <- opt$maximum
        names(lambda) <- "lambda"
        LL <- opt$objective
# get GLS coefficients
        fit <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, I=I,
            weights=Sweights, family=family, out=TRUE, tol.solve=tol.solve)
# create residuals and fitted values (Cressie 1993, p. 564)
	fit$signal_trend <- drop(X %*% fit$coefficients)
	fit$signal_stochastic <- drop(lambda * W %*% (Y - fit$signal_trend))
	fit$fitted.values <- fit$signal_trend + fit$signal_stochastic
	fit$residuals <- drop(Y - fit$fitted.values)
# get null LL
        LL0 <- .opt.fit.spam(lambda=as.numeric(0), Y=Y, X=X, n=n, W=W,
            W_J=W_J, I=I, weights=Sweights, sum_lw=sum_lw, family=family,
            verbose=verbose, tol.solve=tol.solve)
#        weights <- diag(Sweights)
    } else stop("unknown method")
    res <- list(fit=fit, lambda=lambda, LL=LL, LL0=LL0, call=match.call(),
        parameters=(ncol(X)+2), aliased=aliased, method=method,
        zero.policy=zero.policy, weights=weights)
    if (!is.null(na.act))
	res$na.action <- na.act
    if (is.null(llprof)) res$llprof <- llprof
    else {
        res$llprof <- list(lambda=llprof, ll=ll_prof)
    }
    if (zero.policy) {
        zero.regs <- attr(listw$neighbours, 
	    "region.id")[which(card(listw$neighbours) == 0)]
	if (length(zero.regs) > 0)
	    attr(res, "zero.regs") <- zero.regs
	}

    class(res) <- "spautolm"
    res
}

.opt.fit.full <- function(lambda, Y, X, n, W, eig, I, weights, sum_lw,
    family="SAR", verbose=TRUE, tol.solve=.Machine$double.eps) {
# fitting function called from optimize()
#    IlW <- diag(n) - lambda * W
    SSE <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, weights=weights,
        I=I, family=family, out=FALSE, tol.solve=tol.solve)
    s2 <- SSE/n
    if (family == "SMA") {
        if (is.complex(eig)) detIlW <- Re(prod(1/(1 + lambda * eig))) 
        else detIlW <- prod(1/(1 + lambda * eig))
        ldet <- log(detIlW)
    } else {
        if (is.complex(eig)) detIlW <- Re(prod(1 - lambda*eig)) 
        else detIlW <- prod(1 - lambda*eig)
        ldet <- (1/ifelse((length(grep("CAR", family)) != 0), 2, 1)) * 
            log(detIlW)
    }
    ret <- (ldet + (1/2)*sum_lw - ((n/2)*log(2*pi)) - (n/2)*log(s2) - 
        (1/(2*(s2)))*SSE)
    if (verbose)  cat("lambda:", lambda, "function:", ret, "Jacobian", log(detIlW), "SSE", SSE, "\n")
    ret
}

#.opt.fit.SparseM <- function(lambda, Y, X, n, W, W_J, I, weights, sum_lw,
#    family="SAR", verbose=TRUE, cholAlloc, tol.solve=.Machine$double.eps) {
## fitting function called from optimize()
#    SSE <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, weights=weights,
#        I=I, family=family, out=FALSE, tol.solve=tol.solve)
#    s2 <- SSE/n
#    Det <- get("det", "package:SparseM")
#    Jacobian <- log(Det(chol((I - lambda * W_J), nsubmax=cholAlloc$nsubmax, 
#	nnzlmax=cholAlloc$nnzlmax, tmpmax=cholAlloc$tmpmax))^2)
#    gc(FALSE)
#    ret <- ((1/ifelse((length(grep("CAR", family)) != 0), 2, 1))*Jacobian +
#	(1/2)*sum_lw - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
#    if (verbose)  cat("lambda:", lambda, "function", ret, "Jacobian", 
#        Jacobian, "SSE", SSE, "\n")
#    ret
#}

.opt.fit.Matrix <- function(lambda, Y, X, n, W, W_J, I, weights, sum_lw,
    family="SAR", verbose=TRUE, tol.solve=.Machine$double.eps) {
# fitting function called from optimize()
    SSE <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, weights=weights,
        I=I, family=family, out=FALSE, tol.solve=tol.solve)
    s2 <- SSE/n
    CHOL <- try(chol(as((I - lambda * W_J), "dsCMatrix")), silent=TRUE)
    if (class(CHOL) == "try-error") {
        Jacobian <- NA
    } else {
        Jacobian <- sum(2*log(diag(CHOL)))
    }
#    gc(FALSE)
    ret <- ((1/ifelse((length(grep("CAR", family)) != 0), 2, 1))*Jacobian +
	(1/2)*sum_lw - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (verbose)  cat("lambda:", lambda, "function", ret, "Jacobian", Jacobian, "SSE", SSE, "\n")
    ret
}

.opt.fit.spam <- function(lambda, Y, X, n, W, W_J, I, weights, sum_lw,
    family="SAR", verbose=TRUE, tol.solve=.Machine$double.eps) {
# fitting function called from optimize()
    SSE <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, weights=weights,
        I=I, family=family, out=FALSE, tol.solve=tol.solve)
    s2 <- SSE/n
    J1 <- try(determinant((I - lambda * W_J), logarithm=TRUE)$modulus,
        silent=TRUE)
    if (class(J1) == "try-error") {
        Jacobian <- NA
    } else {
        Jacobian <- J1
    }
#    gc(FALSE)
    ret <- ((1/ifelse((length(grep("CAR", family)) != 0), 2, 1))*Jacobian +
	(1/2)*sum_lw - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (verbose)  cat("lambda:", lambda, "function", ret, "Jacobian", Jacobian, "SSE", SSE, "\n")
    ret
}

.SPAR.fit <- function(lambda, Y, X, n, W, weights, I, family,
    out=FALSE, tol.solve=.Machine$double.eps) {
    dmmf <- eval(parse(text=family))
    if (family == "SMA") IlW <- dmmf((I + lambda * W), weights)
    else IlW <- dmmf((I - lambda * W), weights)
    imat <- base:::solve(t(X) %*% as.matrix(IlW %*% X), tol=tol.solve)
    coef <- crossprod(imat, t(X) %*% as.matrix(IlW %*% Y))
    fitted <- X %*% coef
    residuals <- Y - fitted
    SSE <- c(t(residuals) %*% as.matrix(IlW %*% residuals))
    if (!out) return(SSE)

    s2 <- SSE/n
#    var <- s2 * diag(imat)
    coef <- c(coef)
    names(coef) <- colnames(X)
    res <- list(coefficients=coef, SSE=c(SSE), s2=c(s2), imat=imat,
        N=length(residuals))
    res
}

# Simultaneous autoregressive
SAR <- function(IlW, weights) {
    t(IlW) %*% weights %*% IlW
}

# Conditional  autoregressive
CAR <- function(IlW, weights) {
    IlW %*% weights
}

# Spatial moving average
SMA <- function(IlW, weights) {
    IlW <- solve(IlW)
    t(IlW) %*% weights %*% IlW
}


print.spautolm <- function(x, ...) {
	cat("\nCall:\n")
	print(x$call)
	cat("\nCoefficients:\n")
	print(coef(x))
	cat("\nLog likelihood:", logLik(x), "\n")
	invisible(x)
    
}

residuals.spautolm <- function(object, ...) {
	if (is.null(object$na.action))
		object$fit$residuals
	else napredict(object$na.action, object$residuals)
}

fitted.spautolm <- function(object, ...) {
	if (is.null(object$na.action))
		object$fit$fitted.values
	else napredict(object$na.action, object$fitted.values)
}

deviance.spautolm <- function(object, ...) {
	object$SSE
}

coef.spautolm <- function(object, ...) {
	c(object$fit$coefficients, object$lambda)
}


logLik.spautolm <- function(object, ...) {
	LL <- c(object$LL)
	class(LL) <- "logLik"
	N <- object$fit$N
	attr(LL, "nall") <- N
	attr(LL, "nobs") <- N
	attr(LL, "df") <- object$parameters
	LL
}

LR1.spautolm <- function(object)
{
	if (!inherits(object, "spautolm")) stop("Not a spautolm object")
	LLx <- logLik(object)
	LLy <- object$LL0
	statistic <- 2*(LLx - LLy)
	attr(statistic, "names") <- "Likelihood ratio"
	parameter <- 1
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	estimate <- c(LLx, LLy)
	attr(estimate, "names") <- c(paste("Log likelihood of spatial regression fit"), paste("Log likelihood of OLS fit",
		deparse(substitute(y))))
	method <- "Likelihood Ratio diagnostics for spatial dependence"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res
}

summary.spautolm <- function(object, correlation = FALSE, adj.se=FALSE, ...)
{
	N <- object$fit$N
	adj <- ifelse (adj.se, N/(N-length(object$fit$coefficients)), 1) 
	object$fit$s2 <- object$fit$s2*adj
	object$resvar <- object$fit$s2*object$fit$imat
	rownames(object$resvar) <- colnames(object$resvar) <- 
		names(object$fit$coefficients)
	object$adj.se <- adj.se

	object$rest.se <- sqrt(diag(object$resvar))
	object$Coef <- cbind(object$fit$coefficients, object$rest.se, 
		object$fit$coefficients/object$rest.se,
		2*(1-pnorm(abs(object$fit$coefficients/object$rest.se))))
	colnames(object$Coef) <- c("Estimate", "Std. Error", 
		ifelse(adj.se, "t value", "z value"), "Pr(>|z|)")
	if (correlation) {
		object$correlation <- diag((diag(object$resvar))
			^(-1/2)) %*% object$resvar %*% 
			diag((diag(object$resvar))^(-1/2))
		dimnames(object$correlation) <- dimnames(object$resvar)
	}
	object$LR1 <- LR1.spautolm(object)
	rownames(object$Coef) <- names(object$fit$coefficients)
	structure(object, class=c("summary.spautolm", class(object)))
}

print.summary.spautolm <- function(x, digits = max(5, .Options$digits - 3),
	signif.stars = FALSE, ...)
{
	cat("\nCall: ", deparse(x$call),	sep = "", fill=TRUE)
	cat("\nResiduals:\n")
	resid <- residuals(x)
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2) 
		structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
			dimnames(resid)[[2]]))
	else structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
	if (x$zero.policy) {
		zero.regs <- attr(x, "zero.regs")
		if (!is.null(zero.regs))
			cat("\nRegions with no neighbours included:\n",
			zero.regs, "\n")
	}
	cat("\nCoefficients:", x$coeftitle, "\n")
	coefs <- x$Coef
	if (!is.null(aliased <- x$aliased) && any(x$aliased)){
		cat("    (", table(aliased)["TRUE"], 
			" not defined because of singularities)\n", sep = "")
		cn <- names(aliased)
		coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                	colnames(x$Coef)))
            	coefs[!aliased, ] <- x$Coef
	}
	printCoefmat(coefs, signif.stars=signif.stars, digits=digits,
		na.print="NA")
	res <- x$LR1
	cat("\nLambda:", format(signif(x$lambda, digits)),
		"LR test value:", format(signif(res$statistic, digits)),
		"p-value:", format.pval(res$p.value, digits), 
		"\n")
	cat("\nLog likelihood:", logLik(x), "\n")
	if (x$adj.se) cat("Residual variance (sigma squared): ") 
	else cat("ML residual variance (sigma squared): ") 
	cat(format(signif(x$fit$s2, digits)), ", (sigma: ", 
		format(signif(sqrt(x$fit$s2), digits)), ")\n", sep="")
	cat("Number of observations:", x$fit$N, "\n")
	cat("Number of parameters estimated:", x$parameters, "\n")
	cat("AIC: ", format(signif(AIC(x), digits)), "\n", sep="")
    	correl <- x$correlation
    	if (!is.null(correl)) {
        	p <- NCOL(correl)
        	if (p > 1) {
            		cat("\nCorrelation of Coefficients:\n")
                	correl <- format(round(correl, 2), nsmall = 2, 
                  	digits = digits)
                	correl[!lower.tri(correl)] <- ""
                	print(correl[-1, -p, drop = FALSE], quote = FALSE)
            	}
    	}
    	cat("\n")
        invisible(x)
}

