# Copyright 2002-4 by Roger Bivand and Michael Tiefelsdorf,
# with contributions by Danlin Yu
#

localmoran.sad <- function (model, select, nb, glist = NULL, style = "W",
    zero.policy = FALSE, alternative = "greater", spChk=NULL,
    save.Vi = FALSE, tol = .Machine$double.eps^0.5,
    maxiter = 1000, tol.bounds=0.0001) {
# need to impose check on weights TODO!!
    if (class(nb) != "nb") 
        stop(paste(deparse(substitute(nb)), "not an nb object"))
#    if (class(model) != "lm") 
#        stop(paste(deparse(substitute(model)), "not an lm object"))
    clobj <- class(model)
    cond.sad <- FALSE
    n <- length(nb)
    if (clobj == "sarlm") {
    	errorsarLambda <- c(model$lambda)
	if (is.null(errorsarLambda)) stop(paste(deparse(substitute(model)), 
	    "not an error sarlm object"))
	model <- model$lm.model
	W <- nb2mat(nb, glist=glist, style=style, zero.policy=zero.policy)
	sqrtOmega <- solve(diag(n) - errorsarLambda * W)
	cond.sad <- TRUE
     } else if (clobj != "lm")
     	stop(paste(deparse(substitute(model)), "not an lm object"))
    u <- residuals(model)
    if (n != length(u)) 
        stop("objects of different length")
    if (is.null(spChk)) spChk <- get.spChkOption()
    if (spChk && !chkIDs(u, nb2listw(nb, zero.policy=zero.policy)))
	stop("Check of data and weights ID integrity failed")
    if (!(alternative %in% c("greater", "less", "two.sided")))
	stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
    u <- as.vector(u)
    select <- unique(as.integer(select))
    if (any(select < 1 || select > n))
        stop("select out of range")
    utu <- c(t(u) %*% u)
    p <- model$rank
    p1 <- 1:p
    m <- n - p - 2
    XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    X <- model.matrix(terms(model), model.frame(model))
    B <- listw2U(nb2listw(nb, glist=glist, style="B",
	zero.policy=zero.policy))
    D <- NULL
    a <- NULL
    if (style == "W") {
        D <- 1/sapply(B$weights, sum)
    } else if (style == "S") {
        D <- 1 / sqrt(sapply(B$weights, function(x) sum(x^2)))
#        a <- sum(unlist(B$weights))
# correction by Danlin Yu, 25 March 2004
	a <- sum(sapply(B$weights, function(x) sqrt(sum(x^2))))
    } else if (style == "C") a <- sum(unlist(B$weights))
    res <- vector(mode="list", length=length(select))
    for (i in 1:length(select)) {
        Vi <- listw2star(B, select[i], style=style, n, D, a,
	    zero.policy=zero.policy)
        Viu <- lag.listw(Vi, u, zero.policy=TRUE)
	Ii <- c((t(u) %*% Viu) / utu)
	if (cond.sad) {
	    M <- diag(n) - X %*% XtXinv %*% t(X)
    	    ViI <- listw2mat(Vi) - Ii * diag(n)
	    innerTerm <- t(sqrtOmega) %*% M %*% ViI %*% M %*% sqrtOmega
	    evalue <- eigen(innerTerm, only.values=TRUE)$values
#	    is.na(evalue) <- abs(evalue) < 1.0e-7
#	    evalue <- na.omit(evalue)
#	    saveNaAction <- attr(evalue, "na.action")
	    saveNaAction <- NA
	    tau <- c(evalue)
	    e1 <- tau[1]
	    en <- tau[length(tau)]
#	    taumi <- tau
#    	    taumi <- tau - Ii
    	    low <- (1 / (2*tau[length(tau)])) + tol.bounds #+ 0.01
    	    high <- (1 / (2*tau[1])) - tol.bounds #- 0.01
    	    f <- function(omega, tau) {sum(tau/(1 - (2*omega*tau)))}
    	    root <- uniroot(f, lower=low, upper=high, tol=tol, maxiter=maxiter,
        	tau=tau)
    	    omega <- root$root
# 0 should be expectation - maybe use try()
	    if (omega < 0 ) sad.r <- try(-sqrt(sum(log(1 - 2*omega*tau))))
    	    else sad.r <- try(sqrt(sum(log(1 - 2*omega*tau))))
	    if (inherits(sad.r, "try.error")) {
	    	warning (paste("In zone:", select[i], "sad.r not a number"))
	    	sad.r <- sad.u <- sad.p <- NaN
	    } else { 
		sad.u <- omega * sqrt(2*sum(tau^2 / (1 - (2*omega*tau))^2))
    	    	sad.p <- sad.r - ((1/sad.r)*log(sad.r/sad.u))
	    }
	} else {
	    ViX <- lag.listw(Vi, X, zero.policy=TRUE)
	    MViM <- t(X) %*% ViX %*% XtXinv
	    t1 <- -sum(diag(MViM))
	    sumsq.Vi <- function(x) {
            	if (is.null(x)) NA
	    	else sum(x^2)
	    }
	    trVi2 <- sum(sapply(Vi$weights, sumsq.Vi), na.rm=TRUE)
	    t2a <- sum(diag(t(ViX) %*% ViX %*% XtXinv))
	    t2b <- sum(diag(MViM %*% MViM))
	    t2 <- trVi2 - 2*t2a + t2b
	    e1 <- 0.5 * (t1 + sqrt(2*t2 - t1^2))
	    en <- 0.5 * (t1 - sqrt(2*t2 - t1^2))
	    l <- en
	    h <- e1
	    mi <- Ii
	    aroot= m*mi*(l+h-2*mi)+mi*(3*l+3*h-4*mi)-2*l*h
            broot= (m+2)*mi*(l-mi)*(h-mi)
            c1root= l**2 * mi**2 * (m+1)**2 + h**2 * mi**2 * (m+1)**2
            c2root= 2*l*h * (2*l*h - 2*l*mi - 2*h*mi - 2*m*mi**2 -
	    	m**2 * mi**2 + mi**2)
            omega= 0.25*((aroot-sqrt(c1root+c2root))/broot)
	    if (is.nan(omega)) {
	    	warning (paste("In zone:", select[i], "omega not a number"))
	    	sad.r <- sad.u <- sad.p <- NaN
	    } else { 
            	tau <- c(c(e1), rep(0, m), c(en))
	    	taumi <- tau - Ii
            	if (omega < 0 ) sad.r <- -sqrt(sum(log(1 - 2*omega*taumi)))
            	else sad.r <- sqrt(sum(log(1 - 2*omega*taumi)))
            	sad.u <- omega * sqrt(2*sum(taumi^2 / (1 - (2*omega*taumi))^2))
            	sad.p <- sad.r - ((1/sad.r)*log(sad.r/sad.u))
	    }
	}
        if (alternative == "two.sided") p.sad <- 2 * pnorm(abs(sad.p), 
	    lower.tail=FALSE)
        else if (alternative == "greater")
            p.sad <- pnorm(sad.p, lower.tail=FALSE)
        else p.sad <- pnorm(sad.p)
        statistic <- sad.p
        attr(statistic, "names") <- "Saddlepoint approximation"
        p.value <- p.sad
        estimate <- c(Ii)
        attr(estimate, "names") <- "Observed Moran's Ii"
        internal1 <- c(omega, sad.r, sad.u)
        attr(internal1, "names") <- c("omega", "sad.r", "sad.u")
        method <- paste("Saddlepoint approximation for local Moran's I",
            "(Barndorff-Nielsen formula)")
        data.name <- paste("region:", select[i],
	    attr(nb, "region.id")[select[i]],
	    "\n", paste(strwrap(paste("model: ", gsub(" *", " ", 
	    paste(deparse(model$call), sep="", collapse="")))),
	    collapse="\n"),
            "\nneighbours:", deparse(substitute(nb)),
	    "style:", style, "\n")
        obj <- list(statistic = statistic, p.value = p.value,
            estimate = estimate, method = method,
	    alternative = alternative, data.name = data.name,
	    internal1 = internal1, df = (n-p), tau = c(c(e1), c(en)),
	    i = paste(select[i], attr(nb, "region.id")[select[i]]),
#	    if (save.Vi) {Vi = Vi}
	    Vi = if(save.Vi) Vi else NULL)
        class(obj) <- "moransad"
	res[[i]] <- obj
    }
    class(res) <- "localmoransad"
    invisible(res)
}

print.localmoransad <- function(x, ...) {
    extract <- function(x, i) {x[[i]]}
    regnames <- sapply(x, extract, 10)
    est <- sapply(x, extract, 3)
    sad <- sapply(x, extract, 1)
    pval <- sapply(x, extract, 2)
    res <- as.matrix(cbind(est, sad, pval))
    rownames(res) <- regnames
    colnames(res) <- c("Local Morans I", "Saddlepoint", "Pr. (Sad)")
    print(res, ...)
    invisible(res)
}
as.data.frame.localmoransad <- function(x, row.names=NULL, optional=FALSE) {
    n <- length(x)
    res <- matrix(0, nrow=n, ncol=14)
    regnames <- NULL
    if (!is.null(row.names)) 
	if (length(row.names) == n) regnames <- row.names
    if (is.null(regnames))for (i in 1:n) regnames <- c(regnames, x[[i]]$i)
    for (i in 1:n) {
        tau <- x[[i]]$tau
	df <- x[[i]]$df
        tau <- c(tau[1], rep(0, df-2), tau[2])
        max.I <- tau[1]
        min.I <- tau[length(tau)]
        E.I <- sum(tau)/df
        tau <- tau - E.I
        V.I <- (2*sum(tau^2)) / (df*(df+2))
        Z.I <- (x[[i]]$estimate - E.I) / sqrt(V.I)
	if (x[[i]]$alternative == "two.sided") 
	    P.I <- 2 * (1 - pnorm(Z.I))
        else if (x[[i]]$alternative == "greater")
            P.I <- pnorm(Z.I, lower.tail=FALSE)
        else P.I <- pnorm(Z.I)
        Sk.I <- ((8*sum(tau^3))/(df*(df+2)*(df+4))) / (V.I^(3/2))
        Kur.I <- ((48*sum(tau^4) + 12*(sum(tau^2))^2) /
            (df*(df+2)*(df+4)*(df+6))) / (V.I^2)
	res[i,] <- c(x[[i]]$estimate, Z.I, P.I, x[[i]]$statistic,
	    x[[i]]$p.value, E.I, V.I, Sk.I, Kur.I, min.I, max.I,
	    x[[i]]$internal1)
    }
    colnames(res) <- c("Local Morans I", "Stand. dev. (N)", "Pr. (N)",
        "Saddlepoint", "Pr. (Sad)", "Expectation", "Variance",
        "Skewness", "Kurtosis", "Minimum", "Maximum",
        "omega", "sad.r", "sad.u")
    rownames(res) <- regnames
    res <- as.data.frame(res)
    res
}


summary.localmoransad <- function(object, ...) {
    res <- as.data.frame(object)
    class(res) <- c("summary.localmoransad", class(res)) 
    res
}

print.summary.localmoransad <- function(x, ...) {
	print(as.data.frame(x), ...)
	invisible(x)
}



