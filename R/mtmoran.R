# Copyright 2002 by Roger Bivand and Michael Tiefelsdorf
#

lm.morantest.sad <- function (model, listw, zero.policy = FALSE, 
    alternative = "greater", spChk=NULL, tol = .Machine$double.eps^0.5,
    maxiter = 1000) 
{
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (class(model) != "lm") 
        stop(paste(deparse(substitute(model)), "not an lm object"))
    N <- length(listw$neighbours)
    u <- residuals(model)
    if (N != length(u)) 
        stop("objects of different length")
    if (is.null(spChk)) spChk <- get.spChkOption()
    if (spChk && !chkIDs(u, listw))
	stop("Check of data and weights ID integrity failed")
    if (!(alternative %in% c("greater", "less", "two.sided")))
	stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
    u <- as.vector(u)
    listw.U <- listw2U(listw)
    S0 <- sum(unlist(listw.U$weights))
    lu <- lag.listw(listw.U, u, zero.policy = zero.policy)
    Nnn <- N
    if (zero.policy) Nnn <- length(which(card(listw$neighbours) > 0))
    I <- (Nnn/S0) * ((t(u) %*% lu)/(t(u) %*% u))
    p <- model$rank
    p1 <- 1:p
    XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    X <- model.matrix(terms(model), model.frame(model))
    M <- diag(N) - X %*% XtXinv %*% t(X)
    U <- listw2mat(listw.U)
    MVM <- M %*% U %*% M
    MVM <- 0.5 * (t(MVM) + MVM)
    evalue <- eigen(MVM, only.values=TRUE)$values
    idxpos <- (which(abs(evalue) < 1.0e-7)[1]) - 1
    tau <- evalue[1:idxpos]
    tau <- c(tau, evalue[(idxpos+1+p):N])
    taumi <- tau - I
    low <- (1 / (2*taumi[length(taumi)])) + 0.01
    high <- (1 / (2*taumi[1])) - 0.01
    f <- function(omega, taumi) {sum(taumi/(1 - (2*omega*taumi)))}
    root <- uniroot(f, lower=low, upper=high, tol=tol, maxiter=maxiter,
        taumi=taumi)
    omega <- root$root
    if (omega < 0 ) sad.r <- -sqrt(sum(log(1 - 2*omega*taumi)))
    else sad.r <- sqrt(sum(log(1 - 2*omega*taumi)))
    sad.u <- omega * sqrt(2*sum(taumi^2 / (1 - (2*omega*taumi))^2))
    sad.p <- sad.r - ((1/sad.r)*log(sad.r/sad.u))
    if (alternative == "two.sided") p.sad <- 2 * pnorm(abs(sad.p), 
	lower.tail=FALSE)
    else if (alternative == "greater")
        p.sad <- pnorm(sad.p, lower.tail=FALSE)
    else p.sad <- pnorm(sad.p)
    statistic <- sad.p
    attr(statistic, "names") <- "Saddlepoint approximation"
    p.value <- p.sad
    estimate <- c(I)
    attr(estimate, "names") <- "Observed Moran's I"
    internal1 <- c(omega, sad.r, sad.u)
    attr(internal1, "names") <- c("omega", "sad.r", "sad.u")
    internal2 <- unlist(root)[2:4]
    attr(internal2, "names") <- c("f.root", "iter", "estim.prec")
    method <- paste("Saddlepoint approximation for global Moran's I",
        "(Barndorff-Nielsen formula)")
    data.name <- paste("\nmodel:", paste(strwrap(gsub(" *", " ", 
	    paste(deparse(model$call), sep="", collapse=""))), collapse="\n"),
    	    "\nweights: ", deparse(substitute(listw)), "\n", sep="")
    res <- list(statistic = statistic, p.value = p.value,
        estimate = estimate, method = method,
	alternative = alternative, data.name = data.name,
	internal1 = internal1, internal2 = internal2,
	df = (N-p), tau = tau)
    class(res) <- "moransad"
    return(res)
}

print.moransad <- function(x, ...) {
    print.htest(x, ...)
    invisible(x)
}

summary.moransad <- function(object, ...) {
    res <- object
    tau <- object$tau
    df <- object$df
    if (length(tau) == 2) tau <- c(tau[1], rep(0, df-2), tau[2])
    max.I <- tau[1]
    min.I <- tau[length(tau)]
    E.I <- sum(tau)/df
    tau <- tau - E.I
    V.I <- (2*sum(tau^2)) / (df*(df+2))
    Z.I <- (object$estimate - E.I) / sqrt(V.I)
    Sk.I <- ((8*sum(tau^3))/(df*(df+2)*(df+4))) / (V.I^(3/2))
    Kur.I <- ((48*sum(tau^4) + 12*(sum(tau^2))^2) /
        (df*(df+2)*(df+4)*(df+6))) / (V.I^2)
    res$xtra <- c(E.I, V.I, Z.I, Sk.I, Kur.I, min.I, max.I)
    attr(res$xtra, "names") <- c("Expectation", "Variance",
        "Std. deviate", "Skewness", "Kurtosis", "Minimum", "Maximum")
    class(res) <- c("summary.moransad", "moransad")
    res
}

print.summary.moransad <- function(x, ...) {
    print.htest(x, ...)
    print(c(x$xtra, x$internal1), ...)
    if (!is.null(x$internal2)) print(x$internal2, ...)
    invisible(x)
}

