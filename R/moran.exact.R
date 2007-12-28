# Copyright (c) 2007 Markus Reder and Roger Bivand

lm.morantest.exact <- function(model, listw, zero.policy = FALSE, 
    alternative = "greater", spChk=NULL, resfun=weighted.residuals, 
    zero.tol=1.0e-7, Omega=NULL, save.M=NULL, save.U=NULL) 
{
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (class(model) != "lm") 
        stop(paste(deparse(substitute(model)), "not an lm object"))
    N <- length(listw$neighbours)
    u <- resfun(model)
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
    I <- (Nnn/S0) * (crossprod(u, lu) / crossprod(u))
    I_save <- I
    if (!isTRUE(all.equal((Nnn/S0), 1))) I <- I * (S0/Nnn)
    p <- model$rank
    p1 <- 1:p
    nacoefs <- which(is.na(coefficients(model)))
    XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
    X <- model.matrix(terms(model), model.frame(model))
# fixed after looking at TOWN dummy in Boston data
    if (length(nacoefs > 0)) X <- X[,-nacoefs]
    if (!is.null(wts <- weights(model))) {
	X <- sqrt(diag(wts)) %*% X
    }
    M <- diag(N) - X %*% tcrossprod(XtXinv, X)
    U <- listw2mat(listw.U)
    if (is.null(Omega)) {
        MVM <- M %*% U %*% M
        MVM <- 0.5 * (t(MVM) + MVM)
        evalue <- eigen(MVM, only.values=TRUE)$values
        idxpos <- which(abs(evalue) < zero.tol)
        if (length(idxpos) != p)
            warning("number of zero eigenvalues greater than number of variables")
        idxpos <- idxpos[1] - 1
        if (idxpos < 1) stop("invalid first zero eigenvalue index")
        gamma <- evalue[1:idxpos]
        gamma <- c(gamma, evalue[(idxpos+1+p):N])
        res <- exactMoran(I, gamma, alternative=alternative, type="Global")
    } else {
        if (dim(Omega)[1] != N) stop("Omega of different size than data")
        Omega <- chol(Omega)
        res <- exactMoranAlt(I, M, U, Omega, N, alternative=alternative,
            type="Alternative")
    }

    data.name <- paste("\nmodel:", paste(strwrap(gsub("[[:space:]]+", " ", 
	    paste(deparse(model$call), sep="", collapse=""))), collapse="\n"),
    	    "\nweights: ", deparse(substitute(listw)), "\n", sep="")
    res$estimate <- c(I_save)
    res$data.name <- data.name
    res$df <- (N-p)
    if (!is.null(save.M)) res$M <- M
    if (!is.null(save.U)) res$U <- U
    return(res)
}



exactMoran <- function(I, gamma, alternative="greater", type="Global", np2=NULL) {
    if (!(alternative %in% c("greater", "less", "two.sided")))
	stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
    if (type == "Global") integrand <- function(x) {
        sin(0.5 * colSums(atan((gamma-I) %*% t(x)))) /
        (x * apply((1 + (gamma-I)^2 %*% t(x^2))^(1/4), 2, prod))}
    else if (type == "Alternative") integrand <- function(x) {
        sin(0.5 * colSums(atan((gamma) %*% t(x)))) /
        (x * apply((1+(gamma)^2 %*% t(x^2))^(1/4), 2, prod))}
    else if (type == "Local") {
        if (is.null(np2)) stop("Local requires np2")
        min <- gamma[1]
        max <- gamma[2]
        integrand <- function(x) {
            sin(0.5*(atan((min-I)*x)+(np2)*atan((-I)*x) + 
            atan((max-I)*x)))/(x*((1+(min-I)^2*x^2) *
            (1+I^2*x^2)^(np2) * (1+(max-I)^2*x^2))^(1/4))} 
    }
    II <- integrate(integrand, lower=0, upper=Inf)$value
    sd.ex <- qnorm(0.5-II/pi)
    if (alternative == "two.sided") p.v <- 2 * pnorm(sd.ex, 
	lower.tail=FALSE)
    else if (alternative == "greater")
        p.v <- pnorm(sd.ex, lower.tail=FALSE)
    else p.v <- pnorm(sd.ex)
    if (is.nan(p.v) || p.v < 0 || p.v > 1) 
	warning("Out-of-range p-value: reconsider test arguments")
    statistic <- sd.ex
    attr(statistic, "names") <- "Exact standard deviate"
    p.value <- p.v
    estimate <- c(I)
    attr(estimate, "names") <- "Observed Moran's I"
    method <- paste(type, "Moran's I statistic with exact p-value")

    res <- list(statistic = statistic, p.value = p.value,
        estimate = estimate, method = method,
	alternative = alternative, gamma=gamma)
    class(res) <- "moranex"
    return(res)
}

exactMoranAlt <- function(I, M, U, Omega, n, alternative="greater",
    type="Alternative") {
    A <- Omega %*% M %*% (U- diag(I, n)) %*% M %*% t(Omega)
    gamma <- sort(eigen(A)$values)
    obj <- exactMoran(I, gamma, alternative=alternative,
        type=type)
    obj
}


print.moranex <- function(x, ...) {
    class(x) <- c("htest", "moranex")
    print(x, ...)
    invisible(x)
}

