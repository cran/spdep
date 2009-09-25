# Copyright 2005-9 by Roger Bivand
spautolm <- function(formula, data = list(), listw, weights,
    na.action, verbose=FALSE, tol.opt=.Machine$double.eps^(2/3),
    family="SAR", method="full", interval=c(-1,0.999), zero.policy=FALSE,
#    cholAlloc=NULL, 
    super=NULL, Matrix_intern=TRUE, tol.solve=.Machine$double.eps,
    find_interval=FALSE, llprof=NULL) 
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
# bug x for X Bjarke Christensen 090924
	X <- X[,-nacoef]
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
# fix SMA bounds
	full_lower <- eig.range[1]+.Machine$double.eps
	full_upper <- eig.range[2]-.Machine$double.eps
	if (family == "SMA") {
	    full_lower <- -(eig.range[2]-.Machine$double.eps)
	    full_upper <- -(eig.range[1]+.Machine$double.eps)
	}
        I <- diag(n)
# do line search
        dweights <- diag(weights)
        if (!is.null(llprof)) {
            if (length(llprof) == 1)
                llprof <- seq(full_lower,
                    full_upper, length.out=llprof)
            ll_prof <- numeric(length(llprof))
            for (i in seq(along=llprof)) ll_prof[i] <- .opt.fit.full(
                llprof[i], Y=Y, X=X, n=n, W=W, eig=eig, I=I,
                weights=dweights, sum_lw=sum_lw, family=family,
                verbose=verbose, tol.solve=tol.solve)
        }
        opt <- optimize(.opt.fit.full, lower=full_lower,
            upper=full_upper, maximum=TRUE,
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
    }  else if (method == "Matrix") {
        if (family == "SMA") stop("SMA only for full method")
        if (listw$style %in% c("W", "S") && !can.sim)
        stop("Matrix method requires symmetric weights")
        if (listw$style %in% c("B", "C") && 
 	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
	    stop("Matrix method requires symmetric weights")
        if (listw$style == "U") stop("U style not permitted, use C")
        I <- as_dsCMatrix_I(n)
	W <- as_dgRMatrix_listw(listw)
	W <- as(W, "CsparseMatrix")
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
# Jacobian only from symmetric W_J, W can be asymmetric though
        if (listw$style %in% c("W", "S") & can.sim) {
	    W_J <- listw2U_Matrix(similar.listw_Matrix(listw))
#	    similar <- TRUE
	} else W_J <- as_dsTMatrix_listw(listw)
	W_J <- as(W_J, "CsparseMatrix")
	if (!is.null(super)) {
		if (!is.logical(super)) stop("super must be logical")
		Imult <- 2
		if (listw$style == "B") {
                    Imult <- ceiling((2/3)*max(apply(W_J, 1, sum)))
		    interval <- c(-0.5, +0.25)
		} else interval <- c(-2, +1)
                nW_J <- - W_J
		if (super) {
                    super <- !super
                    warning("super=TRUE not yet working")
                }
		pChol <- Cholesky(W_J, super=super, Imult = Imult)
		nChol <- Cholesky(nW_J, super=super, Imult = Imult)
                if (find_interval && Matrix_intern) {
		  ns1 <- last <- 10
		  plambda1 <- seq(sqrt(.Machine$double.eps), interval[2],
                    length.out=ns1)
		
		  while (last >= ns1) {
                   pdet1 <- Matrix:::ldetL2up(nChol, nW_J, 1/plambda1)
		   wp1 <- which(is.finite(pdet1))
		   last <- wp1[length(wp1)]
		   if (last == ns1) plambda1 <- seq(interval[2], 
		       1.5*interval[2], length.out=ns1)
		  }
                  lwp1n <- plambda1[last]
                  lwp2n <- plambda1[last+1]
		  plambda2 <- seq(lwp2n, lwp1n, length.out=ns1)
                  pdet2 <- Matrix:::ldetL2up(nChol, nW_J, 1/plambda2)
		  wp2 <- which(is.finite(pdet2))
                  lwp2n <- plambda2[wp2[length(wp2)]]
		
		  nlambda1 <- seq(interval[1], -sqrt(.Machine$double.eps),
                    length.out=ns1)
		
		  first <- 1
		  while (first == 1) {
                   ndet1 <- Matrix:::ldetL2up(pChol, W_J, 1/(-nlambda1))
		   wn1 <- which(is.finite(ndet1))
		   first <- wn1[1]
		   if (first == 1) plambda1 <- seq(1.5*interval[1], 
			interval[1], length.out=ns1)
		  }

                  lwn1n <- nlambda1[wn1[1]]
                  lwn2n <- nlambda1[wn1[1]-1]
		  nlambda2 <- seq(lwn2n, lwn1n, length.out=ns1)
                  ndet2 <- Matrix:::ldetL2up(pChol, W_J, 1/(-nlambda2))
		  wn2 <- which(is.finite(ndet2))
                  lwn2n <- nlambda2[wn2[1]]
		  interval <- c(lwn2n, lwp2n)
		  if (verbose) cat("using interval:", interval, "\n")
                }
	}
	else nW_J <- pChol <- nChol <- NULL
#	gc(FALSE)
#        Sweights <- as(Diagonal(x=weights), "sparseMatrix")
        Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
	    "CsparseMatrix")
# do line search
        if (!is.null(llprof)) {
            ll_prof <- numeric(length(llprof))
            for (i in seq(along=llprof)) ll_prof[i] <- .opt.fit.Matrix(
                llprof[i], Y=Y, X=X, n=n, W=W, W_J=W_J, I=I,
                weights=Sweights, sum_lw=sum_lw, family=family,
                verbose=verbose, tol.solve=tol.solve, super=super, nW_J=nW_J,
                pChol=pChol, nChol=nChol, Matrix_intern=Matrix_intern)
        }
        opt <- optimize(.opt.fit.Matrix, lower=interval[1],
            upper=interval[2], maximum=TRUE,
            tol = tol.opt, Y=Y, X=X, n=n, W=W, W_J=W_J, I=I,
            weights=Sweights, sum_lw=sum_lw, family=family,
            verbose=verbose, tol.solve=tol.solve, super=super, nW_J=nW_J,
            pChol=pChol, nChol=nChol, Matrix_intern=Matrix_intern)
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
            verbose=verbose, tol.solve=tol.solve, super=super, nW_J=nW_J,
                pChol=pChol, nChol=nChol, Matrix_intern=Matrix_intern)
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
	    W_J <- listw2U_spam(similar.listw_spam(listw))
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


.opt.fit.Matrix <- function(lambda, Y, X, n, W, W_J, I, weights, sum_lw,
    family="SAR", verbose=TRUE, tol.solve=.Machine$double.eps, 
    super=NULL, nW_J=NULL, pChol=NULL, nChol=NULL, Matrix_intern=FALSE) {
# fitting function called from optimize()
    SSE <- .SPAR.fit(lambda=lambda, Y=Y, X=X, n=n, W=W, weights=weights,
        I=I, family=family, out=FALSE, tol.solve=tol.solve)
    s2 <- SSE/n
    .f <- if (package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1
    if (isTRUE(all.equal(lambda, 0))) {
        Jacobian <- lambda
    } else {
        if (is.null(super))
            Jacobian <- determinant(I - lambda * W_J, logarithm=TRUE)$modulus
        else {
            if (lambda > 0) {
                if (Matrix_intern) 
                    detTRY <- try(Matrix:::ldetL2up(nChol, nW_J, 1/lambda),
                        silent=TRUE)
                else
                    detTRY <- try(c(determinant(update(nChol, nW_J, 
                        1/lambda))$modulus), silent=TRUE)
                if (class(detTRY) == "try-error") {
                    Jacobian <- NaN
                } else {
                    Jacobian <- n * log(lambda) + (.f * detTRY)
                }
            } else {
                if (Matrix_intern) 
                    detTRY <- try(Matrix:::ldetL2up(pChol, W_J, 1/(-lambda)),
                        silent=TRUE)
                else
                    detTRY <- try(c(determinant(update(pChol, W_J, 
                        1/(-lambda)))$modulus), silent=TRUE)
                if (class(detTRY) == "try-error") {
                    Jacobian <- NaN
                } else {
                    Jacobian <- n * log(-(lambda)) + (.f * detTRY)
                }
            }
        }
    }
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
	else napredict(object$na.action, object$fit$residuals)
}

fitted.spautolm <- function(object, ...) {
	if (is.null(object$na.action))
		object$fit$fitted.values
	else napredict(object$na.action, object$fit$fitted.values)
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

