# Copyright 2005-2010 by Roger Bivand
spautolm <- function(formula, data = list(), listw, weights,
    na.action, family="SAR", method="full", verbose=NULL,
    interval=NULL, zero.policy=NULL, tol.solve=.Machine$double.eps, llprof=NULL,
    control=list()) {
    timings <- list()
    .ptime_start <- proc.time()
    con <- list(tol.opt=.Machine$double.eps^(2/3), 
       Imult=2, super=NULL, cheb_q=5, MC_p=16, MC_m=30, spamPivot="MMD",
       in_coef=0.1)
    nmsC <- names(con)
    con[(namc <- names(control))] <- control
    if (length(noNms <- namc[!namc %in% nmsC])) 
        warning("unknown names in control: ", paste(noNms, collapse = ", "))
    if (!inherits(listw, "listw")) 
        stop("No neighbourhood list")
    if (is.null(verbose)) verbose <- get("verbose", envir = .spdepOptions)
    stopifnot(is.logical(verbose))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))

    if (family == "SMA" && method != "full") stop("SMA only for full method")
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
    if (NROW(X) != length(listw$neighbours))
	 stop("Input data and neighbourhood list have different dimensions")
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
    can.sim <- FALSE
    if (listw$style %in% c("W", "S")) 
	can.sim <- can.be.simmed(listw)

    sum_lw <- sum(log(weights))
#    env <- new.env(parent=globalenv())
    env <- new.env()
    assign("Y", Y, envir=env)
    assign("X", X, envir=env)
    assign("n", n, envir=env)
    assign("weights", weights, envir=env)
    assign("can.sim", can.sim, envir=env)
    assign("family", family, envir=env)
    assign("method", method, envir=env)
    assign("verbose", verbose, envir=env)
    assign("listw", listw, envir=env)
    assign("sum_lw", sum_lw, envir=env)
    timings[["set_up"]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()


    if (method == "full") {
        eigen_setup(env)
        er <- get("eig.range", envir=env)
        if (is.null(interval)) 
            interval <- c(er[1]+.Machine$double.eps, 
                er[2]-.Machine$double.eps)
# fix SMA bounds
	if (family == "SMA") interval <- -rev(interval)
# spatial weights matrix
        W <- listw2mat(listw)
	attr(W, "dimnames") <- NULL
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
        I <- diag(n)
        Sweights <- diag(weights)
        assign("W", W, envir=env)
        assign("I", I, envir=env)
        assign("Sweights", Sweights, envir=env)
    }  else if (method == "Matrix") {
        if (listw$style %in% c("W", "S") && !can.sim)
        stop("Matrix method requires symmetric weights")
        if (listw$style %in% c("B", "C") && 
 	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
	    stop("Matrix method requires symmetric weights")
        if (listw$style == "U") stop("U style not permitted, use C")
	W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
        assign("W", W, envir=env)
        Imult <- con$Imult
        if (is.null(interval)) {
	    if (listw$style == "B") {
                Imult <- ceiling((2/3)*max(sapply(listw$weights, sum)))
                interval <- c(-0.5, +0.25)
            } else interval <- c(-1, 0.999)
        }
# FIXME
        if (is.null(con$super)) con$super <- as.logical(NA)
        Matrix_setup(env, Imult, con$super)
        I <- as_dsCMatrix_I(n)
        assign("I", I, envir=env)
        Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
	    "CsparseMatrix")
        assign("Sweights", Sweights, envir=env)
    }  else if (method == "Matrix_J") {
        if (listw$style %in% c("W", "S") && !can.sim)
        stop("Matrix method requires symmetric weights")
        if (listw$style %in% c("B", "C") && 
 	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
	    stop("Matrix method requires symmetric weights")
        if (listw$style == "U") stop("U style not permitted, use C")
	W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
        assign("W", W, envir=env)
        Imult <- con$Imult
        if (is.null(interval)) {
	    if (listw$style == "B") {
                Imult <- ceiling((2/3)*max(sapply(listw$weights, sum)))
                interval <- c(-0.5, +0.25)
            } else interval <- c(-1, 0.999)
        }
# FIXME
        if (is.null(con$super)) con$super <- FALSE
        Matrix_J_setup(env, super=con$super)
        I <- as_dsCMatrix_I(n)
        assign("I", I, envir=env)
        Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
	    "CsparseMatrix")
        assign("Sweights", Sweights, envir=env)


    }  else if (method == "spam") {
        if (!require(spam)) stop("spam not available")
        if (listw$style %in% c("W", "S") && !can.sim)
        stop("spam method requires symmetric weights")
        if (listw$style %in% c("B", "C", "U") && 
 	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
	    stop("spam method requires symmetric weights")
	W <- as.spam.listw(listw)
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
        assign("W", W, envir=env)
# Jacobian only from symmetric W_J, W can be asymmetric though
        spam_setup(env, pivot=con$spamPivot)
        Sweights <- diag.spam(x=weights, n, n)
        assign("Sweights", Sweights, envir=env)
        if (is.null(interval)) interval <- c(-1, 0.999)
    }  else if (method == "spam_update") {
        if (!require(spam)) stop("spam not available")
        if (listw$style %in% c("W", "S") && !can.sim)
        stop("spam method requires symmetric weights")
        if (listw$style %in% c("B", "C", "U") && 
 	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
	    stop("spam method requires symmetric weights")
	W <- as.spam.listw(listw)
        if (family == "CAR") if (!isTRUE(all.equal(W, t(W))))
	    warning("Non-symmetric spatial weights in CAR model")
        assign("W", W, envir=env)
# Jacobian only from symmetric W_J, W can be asymmetric though
        spam_update_setup(env, in_coef=con$in_coef, pivot=con$spamPivot)
        Sweights <- diag.spam(x=weights, n, n)
        assign("Sweights", Sweights, envir=env)
        if (is.null(interval)) interval <- c(-1, 0.999)
    } else if (method == "LU") {
        LU_setup(env)
        W <- get("W", envir=env)
        Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
	    "CsparseMatrix")
        assign("Sweights", Sweights, envir=env)
        if (is.null(interval)) interval <- c(-1,0.999)
    } else if (method == "MC") {
	if (!listw$style %in% c("W"))
	    stop("MC method requires row-standardised weights")
        mcdet_setup(env, p=con$MC_p, m=con$MC_m)
        I <- as_dsCMatrix_I(n)
        assign("I", I, envir=env)
        Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
	    "CsparseMatrix")
        assign("Sweights", Sweights, envir=env)
        W <- get("W", envir=env)
        if (is.null(interval)) interval <- c(-1,0.999)
    } else if (method == "Chebyshev") {
	if (listw$style %in% c("W", "S") && !can.sim)
	    stop("Chebyshev method requires symmetric weights")
	if (listw$style %in% c("B", "C", "U") && 
	    !(is.symmetric.glist(listw$neighbours, listw$weights)))
	        stop("Chebyshev method requires symmetric weights")
        cheb_setup(env, q=con$cheb_q)
        I <- as_dsCMatrix_I(n)
        assign("I", I, envir=env)
        Sweights <- as(as(Diagonal(x=weights), "symmetricMatrix"), 
	    "CsparseMatrix")
        assign("Sweights", Sweights, envir=env)
        W <- get("W", envir=env)
        if (is.null(interval)) {
 	    if (listw$style == "B") interval <- c(-0.5, +0.25)
            else interval <- c(-1,0.999)
        }
    } else stop("unknown method")

    nm <- paste(method, "set_up", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

    if (!is.null(llprof)) {
        if (length(llprof) == 1L)
            llprof <- seq(interval[1], interval[2], length.out=llprof)
        ll_prof <- numeric(length(llprof))
        for (i in seq(along=llprof)) 
            ll_prof[i] <- .opt.fit(llprof[i], env=env, tol.solve=tol.solve)
        nm <- paste(method, "profile", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
    }

    opt <- optimize(.opt.fit, interval=interval, maximum=TRUE,
        tol = con$tol.opt, env=env, tol.solve=tol.solve)
    lambda <- opt$maximum
    if (isTRUE(all.equal(lambda, interval[1])) ||
        isTRUE(all.equal(lambda, interval[2]))) 
        warning("lambda on interval bound - results should not be used")
    names(lambda) <- "lambda"
    LL <- opt$objective
    nm <- paste(method, "opt", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    .ptime_start <- proc.time()

# get GLS coefficients
    fit <- .SPAR.fit(lambda=lambda, env, out=TRUE, tol.solve=tol.solve)
# create residuals and fitted values (Cressie 1993, p. 564)
    fit$signal_trend <- drop(X %*% fit$coefficients)
    fit$signal_stochastic <- drop(lambda * W %*% (Y - fit$signal_trend))
    fit$fitted.values <- fit$signal_trend + fit$signal_stochastic
    fit$residuals <- drop(Y - fit$fitted.values)

# get null LL
    LL0 <- .opt.fit(lambda=0, env, tol.solve=tol.solve)
# NK null
    LLNullLlm <- logLik(lm(Y ~ 1, weights=weights))
    nm <- paste(method, "output", sep="_")
    timings[[nm]] <- proc.time() - .ptime_start
    rm(env)
    GC <- gc()
    res <- list(fit=fit, lambda=lambda, LL=LL, LL0=LL0, call=match.call(),
        parameters=(ncol(X)+2), aliased=aliased, method=method,
        zero.policy=zero.policy, weights=weights, interval=interval,
        timings=do.call("rbind", timings)[, c(1, 3)], LLNullLlm=LLNullLlm)
    if (!is.null(na.act))
	res$na.action <- na.act
    if (is.null(llprof)) res$llprof <- llprof
    else {
        res$llprof <- list(lambda=llprof, ll=ll_prof)
    }
    if (zero.policy) {
        zero.regs <- attr(listw$neighbours, 
	    "region.id")[which(card(listw$neighbours) == 0)]
	if (length(zero.regs) > 0L)
	    attr(res, "zero.regs") <- zero.regs
	}

    class(res) <- "spautolm"
    res
}

.opt.fit <- function(lambda, env, tol.solve=.Machine$double.eps) {
# fitting function called from optimize()
    SSE <- .SPAR.fit(lambda=lambda, env=env, out=FALSE, tol.solve=tol.solve)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- do_ldet(lambda, env)
    det <- ifelse(get("family", envir=env) == "CAR", 0.5*ldet, ldet)
    ret <- (det + (1/2)*get("sum_lw", envir=env) - ((n/2)*log(2*pi)) - 
        (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env))  cat("lambda:", lambda, "function:", ret, "Jacobian", ldet, "SSE", SSE, "\n")
    ret
}


.SPAR.fit <- function(lambda, env, out=FALSE, tol.solve=.Machine$double.eps) {
    dmmf <- eval(parse(text=get("family", envir=env)))
    if (get("family", envir=env) == "SMA") IlW <- dmmf((get("I", envir=env) + 
        lambda * get("W", envir=env)), get("Sweights", envir=env))
    else IlW <- dmmf((get("I", envir=env) - lambda * get("W", envir=env)), 
        get("Sweights", envir=env))
    X <- get("X", envir=env)
    Y <- get("Y", envir=env)
    imat <- base:::solve(crossprod(X, as.matrix(IlW %*% X)), tol=tol.solve)
    coef <- crossprod(imat, crossprod(X, as.matrix(IlW %*% Y)))
    fitted <- X %*% coef
    residuals <- Y - fitted
    SSE <- c(crossprod(residuals, as.matrix(IlW %*% residuals)))
    if (!out) return(SSE)

    n <- get("n", envir=env)
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
        if (isTRUE(all.equal(x$lambda, x$interval[1])) ||
            isTRUE(all.equal(x$lambda, x$interval[2]))) 
            warning("lambda on interval bound - results should not be used")
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

summary.spautolm <- function(object, correlation = FALSE, adj.se=FALSE,
 Nagelkerke=FALSE, ...) {
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
        if (Nagelkerke) {
            nk <- NK.sarlm(object)
            if (!is.null(nk)) object$NK <- nk
        }
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
        if (isTRUE(all.equal(x$lambda, x$interval[1])) ||
            isTRUE(all.equal(x$lambda, x$interval[2]))) 
            warning("lambda on interval bound - results should not be used")
	cat("\nResiduals:\n")
	resid <- residuals(x)
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2L) 
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
        if (!is.null(x$NK)) cat("Nagelkerke pseudo-R-squared:",
            format(signif(x$NK, digits)), "\n")
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

