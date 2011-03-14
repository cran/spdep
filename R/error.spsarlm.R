# Copyright 1998-2011 by Roger Bivand (non-W styles Rein Halbersma)
#

errorsarlm <- function(formula, data = list(), listw, na.action, etype="error",
	method="eigen", quiet=NULL, zero.policy=NULL, interval=NULL, 
	tol.solve=1.0e-10, trs=NULL, control=list()) {
        timings <- list()
        .ptime_start <- proc.time()
        con <- list(tol.opt=.Machine$double.eps^0.5, returnHcov=TRUE,
            pWOrder=250, fdHess=NULL, optimHess=FALSE, LAPACK=FALSE,
           compiled_sse=FALSE, Imult=2, cheb_q=5, MC_p=16, MC_m=30,
           super=NULL, spamPivot="MMD", in_coef=0.1, type="MC",
           correct=TRUE, trunc=TRUE)
        nmsC <- names(con)
        con[(namc <- names(control))] <- control
        if (length(noNms <- namc[!namc %in% nmsC])) 
            warning("unknown names in control: ", paste(noNms, collapse = ", "))
        if (is.null(quiet)) quiet <- !get("verbose", env = .spdepOptions)
        stopifnot(is.logical(quiet))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", env = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!inherits(listw, "listw")) stop("No neighbourhood list")
        if (is.null(con$fdHess)) con$fdHess <- method != "eigen"
        stopifnot(is.logical(con$fdHess))
        stopifnot(is.logical(con$optimHess))
        stopifnot(is.logical(con$LAPACK))
#        stopifnot(is.logical(con$super))
        stopifnot(is.logical(con$compiled_sse))
        stopifnot(is.character(con$spamPivot))
	switch(etype, error = if (!quiet)
                cat("\nSpatial autoregressive error model\n"),
	    emixed = if (!quiet)
                cat("\nSpatial mixed autoregressive error model\n"),
	    stop("\nUnknown model type\n"))
	can.sim <- FALSE
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	y <- model.response(mf, "numeric")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")
	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	n <- NROW(x)
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (etype != "error") {
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
# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
        LL_null_lm <- NULL
	if ("(Intercept)" %in% colnames(x)) LL_null_lm <- logLik(lm(y ~ 1))
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
# added no intercept Guillaume Blanchet 091103
	if (m > 1 || (m == 1 && K == 1)) {
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

        env <- new.env(parent=globalenv())
        assign("y", y, envir=env)
        assign("x", x, envir=env)
        assign("wy", wy, envir=env)
        assign("WX", WX, envir=env)
        assign("n", n, envir=env)
        assign("p", m, envir=env)
        assign("verbose", !quiet, envir=env)
        assign("family", "SAR", envir=env)
        assign("compiled_sse", con$compiled_sse, envir=env)
        assign("first_time", TRUE, envir=env)
        assign("LAPACK", con$LAPACK, envir=env)
        assign("can.sim", can.sim, envir=env)
        assign("listw", listw, envir=env)
        assign("similar", FALSE, envir=env)
        timings[["set_up"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()

	if (!quiet) cat(paste("\nJacobian calculated using "))
	switch(method,
		eigen = {
                    if (!quiet) cat("neighbourhood matrix eigenvalues\n")
                    eigen_setup(env)
                    er <- get("eig.range", envir=env)
                    if (is.null(interval)) 
                        interval <- c(er[1]+.Machine$double.eps, 
                            er[2]-.Machine$double.eps)
                },
	        Matrix = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		        stop("Matrix method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		        stop("Matrix method requires symmetric weights")
		    if (!quiet) cat("sparse matrix Cholesky decomposition\n")
	            Imult <- con$Imult
                    if (is.null(interval)) {
	                if (listw$style == "B") {
                            Imult <- ceiling((2/3) * max(sapply(listw$weights,
                                sum)))
	                    interval <- c(-0.5, +0.25)
	                } else interval <- c(-1, 0.999)
                    }
                    if (is.null(con$super)) con$super <- as.logical(NA)
                    Matrix_setup(env, Imult, con$super)
                    W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        	    I <- as_dsCMatrix_I(n)
		},
	        Matrix_J = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		        stop("Matrix method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		        stop("Matrix method requires symmetric weights")
		    if (!quiet) cat("sparse matrix Cholesky decomposition\n")
                    if (is.null(interval)) {
	                if (listw$style == "B") {
	                    interval <- c(-0.5, +0.25)
	                } else interval <- c(-1, 0.999)
                    }
                    if (is.null(con$super)) con$super <- FALSE
                    Matrix_J_setup(env, super=con$super)
                    W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        	    I <- as_dsCMatrix_I(n)
		},
	        spam = {
                    if (!require(spam)) stop("spam not available")
		    if (listw$style %in% c("W", "S") && !can.sim)
		        stop("spam method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		        stop("spam method requires symmetric weights")
		    if (!quiet) cat("sparse matrix Cholesky decomposition\n")
                    spam_setup(env)
                    W <- as.spam.listw(get("listw", envir=env))
                    if (is.null(interval)) interval <- c(-1,0.999)
		},
	        spam_update = {
                    if (!require(spam)) stop("spam not available")
		    if (listw$style %in% c("W", "S") && !can.sim)
		        stop("spam method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		        stop("spam method requires symmetric weights")
		    if (!quiet) cat("sparse matrix Cholesky decomposition\n")
                    spam_update_setup(env, in_coef=con$in_coef,
                        pivot=con$spamPivot)
                    W <- as.spam.listw(get("listw", envir=env))
                    if (is.null(interval)) interval <- c(-1,0.999)
		},
                Chebyshev = {
		    if (listw$style %in% c("W", "S") && !can.sim)
		        stop("Chebyshev method requires symmetric weights")
		    if (listw$style %in% c("B", "C", "U") && 
			!(is.symmetric.glist(listw$neighbours, listw$weights)))
		        stop("Chebyshev method requires symmetric weights")
		    if (!quiet) cat("sparse matrix Chebyshev approximation\n")
                    cheb_setup(env, q=con$cheb_q)
                    W <- get("W", envir=env)
        	    I <- as_dsCMatrix_I(n)
                    if (is.null(interval)) interval <- c(-1,0.999)
                },
                MC = {
		    if (!listw$style %in% c("W"))
		        stop("MC method requires row-standardised weights")
		    if (!quiet) cat("sparse matrix Monte Carlo approximation\n")
                    mcdet_setup(env, p=con$MC_p, m=con$MC_m)
                    W <- get("W", envir=env)
        	    I <- as_dsCMatrix_I(n)
                    if (is.null(interval)) interval <- c(-1,0.999)
                },
                LU = {
		    if (!quiet) cat("sparse matrix LU decomposition\n")
                    LU_setup(env)
                    W <- get("W", envir=env)
                    I <- get("I", envir=env)
                    if (is.null(interval)) interval <- c(-1,0.999)
                },
                moments = {
		    if (!quiet) cat("Smirnov/Anselin (2009)", 
                        "trace approximation\n")
                    moments_setup(env, trs=trs, m=con$MC_m, p=con$MC_p,
                        type=con$type, correct=con$correct, trunc=con$trunc)
                    W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
        	    I <- as_dsCMatrix_I(n)
                    if (is.null(interval)) interval <- c(-1,0.999)
                },
		stop("...\n\nUnknown method\n"))

        nm <- paste(method, "set_up", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        if (con$compiled_sse) {
             ptr <- .Call("opt_error_init", PACKAGE="spdep")
             assign("ptr", ptr, envir=env)
        }
	opt <- optimize(sar.error.f, interval=interval, 
		maximum=TRUE, tol=con$tol.opt, env=env)
	lambda <- opt$maximum
        if (isTRUE(all.equal(lambda, interval[1])) ||
            isTRUE(all.equal(lambda, interval[2]))) 
            warning("lambda on interval bound - results should not be used")
	names(lambda) <- "lambda"
	LL <- opt$objective
        if (con$compiled_sse) {
             .Call("opt_error_free", get("ptr", envir=env), PACKAGE="spdep")
        }
        nm <- paste(method, "opt", sep="_")
        timings[[nm]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
	lm.target <- lm(I(y - lambda*wy) ~ I(x - lambda*WX) - 1)
	r <- as.vector(residuals(lm.target))
	fit <- as.vector(y - r)
	p <- lm.target$rank
	SSE <- deviance(lm.target)
	s2 <- SSE/n
	rest.se <- (summary(lm.target)$coefficients[,2])*sqrt((n-p)/n)
	coef.lambda <- coefficients(lm.target)
	names(coef.lambda) <- xcolnames
	lm.model <- lm(y ~ x - 1)
	ase <- FALSE
	lambda.se <- NULL
	LMtest <- NULL
	asyvar1 <- FALSE
        Hcov <- NULL
        timings[["coefs"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        assign("first_time", TRUE, envir=env)
	if (method == "eigen") {
		tr <- function(A) sum(diag(A))
		W <- listw2mat(listw)
		A <- solve(diag(n) - lambda*W)
		WA <- W %*% A
		asyvar <- matrix(0, nrow=2+p, ncol=2+p)
		asyvar[1,1] <- n / (2*(s2^2))
#		asyvar[1,1] <- n / (2*(s2))
		asyvar[2,1] <- asyvar[1,2] <- tr(WA) / s2
#		asyvar[2,1] <- asyvar[1,2] <- tr(WA)
		asyvar[2,2] <- tr(WA %*% WA) + tr(crossprod(WA))
# bug found 100224 German Muchnik Izon
#		asyvar[3:(p+2),3:(p+2)] <- s2*(t(x - lambda*WX) %*% 
                xl <- (x - lambda*WX)
#		asyvar[3:(p+2),3:(p+2)] <- crossprod(xl)
		asyvar[3:(p+2),3:(p+2)] <- crossprod(xl)/s2
		asyvar1 <- solve(asyvar, tol=tol.solve)
		rownames(asyvar1) <- colnames(asyvar1) <- 
			c("sigma", "lambda", xcolnames)
		
		lambda.se <- sqrt(asyvar1[2,2])
#		lambda.se <- sqrt(s2*asyvar1[2,2])
                timings[["eigen_se"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                if (con$returnHcov) {
                    pp <- lm.model$rank
                    p1 <- 1L:pp
                    R <- chol2inv(lm.model$qr$qr[p1, p1, drop = FALSE])
                    B <- tcrossprod(R, x) %*% A
                    A <- solve(diag(n) - lambda*t(W))
                    C <- A %*% x %*% R
                    Hcov <- B %*% C
                    attr(Hcov, "method") <- method
                    timings[["eigen_hcov"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
		ase <- TRUE
	} else {
                if (con$returnHcov) {
                    pp <- lm.model$rank
                    p1 <- 1L:pp
                    R <- chol2inv(lm.model$qr$qr[p1, p1, drop = FALSE])
                    B <- tcrossprod(R, x)
                    B1 <- as(powerWeights(W=W, rho=lambda, order=con$pWOrder,
                        X=B, tol=tol.solve), "matrix")
                    C <- x %*% R
                    C1 <- as(powerWeights(W=t(W), rho=lambda, order=con$pWOrder,
                        X=C, tol=tol.solve), "matrix")
                    Hcov <- B1 %*% C1
                    attr(Hcov, "method") <- method
                    timings[["sparse_hcov"]] <- proc.time() - .ptime_start
                    .ptime_start <- proc.time()
                }
        }
        if (con$fdHess) {
            coefs <- c(lambda, coef.lambda)
            if (con$compiled_sse) {
                ptr <- .Call("hess_error_init", PACKAGE="spdep")
                assign("ptr", ptr, envir=env)
            }
            fdHess <- getVmate(coefs, env, s2, trs, tol.solve=tol.solve,
                optim=con$optimHess)
            if (con$compiled_sse) {
                .Call("hess_error_free", get("ptr", envir=env),
                    PACKAGE="spdep")
            }
            if (is.null(trs)) {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("lambda", colnames(x))
                if (method != "eigen") {
                    lambda.se <- sqrt(fdHess[1, 1])
                }
            } else {
                rownames(fdHess) <- colnames(fdHess) <- 
                    c("sigma2", "lambda", colnames(x))
                if (method != "eigen") {
                    lambda.se <- sqrt(fdHess[2, 2])
                }
            }
            nm <- paste(method, "fdHess", sep="_")
            timings[[nm]] <- proc.time() - .ptime_start
        }
	call <- match.call()
	names(r) <- names(y)
	names(fit) <- names(y)
	ret <- structure(list(type="error", lambda=lambda,
		coefficients=coef.lambda, rest.se=rest.se, 
		LL=LL, s2=s2, SSE=SSE, parameters=(m+2), lm.model=lm.model, 
		method=method, call=call, residuals=r, lm.target=lm.target,
		opt=opt, fitted.values=fit, ase=ase, formula=formula,
		se.fit=NULL, resvar=asyvar1, similar=get("similar", envir=env),
		lambda.se=lambda.se, LMtest=LMtest, zero.policy=zero.policy, 
		aliased=aliased, LLNullLlm=LL_null_lm, Hcov=Hcov,
                interval=interval, fdHess=fdHess,
                optimHess=con$optimHess, insert=!is.null(trs),
                timings=do.call("rbind", timings)[, c(1, 3)]),
                class=c("sarlm"))
        rm(env)
        GC <- gc()
	if (zero.policy) {
		zero.regs <- attr(listw$neighbours, 
			"region.id")[which(card(listw$neighbours) == 0)]
		if (length(zero.regs) > 0L)
			attr(ret, "zero.regs") <- zero.regs
	}
	if (!is.null(na.act))
		ret$na.action <- na.act
	ret
}

sar_error_sse <- function(lambda, env) {
    if (get("compiled_sse", envir=env)) {
        ft <- get("first_time", envir=env)
        SSE <- .Call("R_ml_sse_env", env, lambda, PACKAGE="spdep")
        if (ft) assign("first_time", FALSE, envir=env)
    } else {
        yl <- get("y", envir=env) - lambda * get("wy", envir=env)
        xl <- get("x", envir=env) - lambda * get("WX", envir=env)
	xl.q <- qr.Q(qr(xl, LAPACK=get("LAPACK", envir=env)))
	xl.q.yl <- crossprod(xl.q, yl)
	SSE <- crossprod(yl) - crossprod(xl.q.yl)
    }
    SSE
}


sar.error.f <- function(lambda, env) {
    SSE <- sar_error_sse(lambda, env)
    n <- get("n", envir=env)
    s2 <- SSE/n
    ldet <- do_ldet(lambda, env)
    ret <- (ldet - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*(s2)))*SSE)
    if (get("verbose", envir=env)) cat("lambda:", lambda, " function:", ret, " Jacobian:", ldet, " SSE:", SSE, "\n")
    ret
}


