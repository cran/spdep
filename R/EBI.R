# Copyright 2002-2003 by Roger Bivand and Marilia Carvalho
#

EBImoran <- function (z, listw, nn, S0, zero.policy = FALSE) 
{
    zm <- mean(z)
    zz <- sum((z - zm)^2)
    lz <- lag.listw(listw, z, zero.policy = zero.policy)
    EBI <- (nn/S0) * ((t(z) %*% lz)/zz)
    res <- EBI
    res
}
EBImoran.mc <- function (n, x, listw, nsim, zero.policy = FALSE,
 alternative = "greater", spChk = NULL) 
{
    if (!inherits(listw, "listw")) 
        stop(paste(deparse(substitute(listw)), "is not a listw object"))
    if (missing(nsim)) 
        stop("nsim must be given")
    m <- length(listw$neighbours)
    if (m != length(x)) 
        stop("objects of different length")
    if (m != length(n)) 
        stop("objects of different length")
    if (is.null(spChk)) 
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(x, listw)) 
        stop("Check of data and weights ID integrity failed")
    if (spChk && !chkIDs(n, listw)) 
        stop("Check of data and weights ID integrity failed")
    if (nsim > gamma(m + 1)) 
        stop("nsim too large for this number of observations")
    S0 <- Szero(listw)
    EB <- EBest(n, x)
    p <- EB$raw
    b <- attr(EB, "parameters")$b
    a <- attr(EB, "parameters")$a
    v <- a + (b/x)
    v[v < 0] <- b/x
    z <- (p - b)/sqrt(v)
    res <- numeric(length = nsim + 1)
    for (i in 1:nsim) res[i] <- EBImoran(sample(z), listw, m, 
        S0, zero.policy)
    res[nsim + 1] <- EBImoran(z, listw, m, S0, zero.policy)
    rankres <- rank(res)
    zrank <- rankres[length(res)]
    diff <- nsim - zrank
    diff <- ifelse(diff > 0, diff, 0)
    pval <- (diff + 1)/(nsim + 1)
    if (alternative == "less") 
        pval <- 1 - pval
    else if (alternative == "two.sided") 
        pval <- 2 * pval
    statistic <- res[nsim + 1]
    names(statistic) <- "statistic"
    parameter <- zrank
    names(parameter) <- "observed rank"
    method <- "Monte-Carlo simulation of Empirical Bayes Index"
    data.name <- paste("cases: ", deparse(substitute(n)), ", risk population: ", 
        deparse(substitute(x)), "\nweights: ", deparse(substitute(listw)), 
        "\nnumber of simulations + 1: ", nsim + 1, "\n", sep = "")
    lres <- list(statistic = statistic, parameter = parameter, 
        p.value = pval, alternative = alternative, method = method, 
        data.name = data.name, res = res, z = z)
    class(lres) <- c("htest", "mc.sim")
    lres
}

probmap <- function(n, x) {
    if (!is.numeric(x)) 
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (!is.numeric(n)) 
        stop(paste(deparse(substitute(n)), "is not a numeric vector"))
    if (any(is.na(x))) 
        stop("NA in at risk population")
    if (any(is.na(n))) 
        stop("NA in cases")
    if (any(x < 0)) 
        stop("negative risk population")
    if (any(n < 0)) 
        stop("negative number of cases")
    p <- n/x
    nsum <- sum(n)
    xsum <- sum(x)
    b <- nsum/xsum
    expCount <- x*b
    relRisk <- 100*(n/expCount)
    pmap <- ppois(n, expCount)
    res <- data.frame(raw=p, expCount=expCount, relRisk=relRisk, pmap=pmap)
    res
}


EBest <- function(n, x) {
    if (!is.numeric(x)) 
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (!is.numeric(n)) 
        stop(paste(deparse(substitute(n)), "is not a numeric vector"))
    if (any(is.na(x))) 
        stop("NA in at risk population")
    if (any(is.na(n))) 
        stop("NA in cases")
    if (any(x < .Machine$double.eps)) 
        stop("non-positive risk population")
    if (any(n < 0)) 
        stop("negative number of cases")
    m <- length(n)
    p <- n/x
    nsum <- sum(n)
    xsum <- sum(x)
    b <- nsum/xsum
    s2 <- sum(x * (((p - b)^2)/xsum))
    a <- s2 - (b/(xsum/m))
    if (a < 0) a <- 0
    est <- b + (a*(p - b)) / (a + (b/x))
    res <- data.frame(raw=p, estmm=est)
    attr(res, "parameters") <- list(a=a, b=b)
    res
}

EBlocal <- function(n, x, nb, zero.policy = FALSE, spChk = NULL) {
    if (class(nb) != "nb") 
        stop(paste(deparse(substitute(nb)), "is not an nb object"))
    m <- length(nb)
    if (m != length(x)) 
        stop("objects of different length")
    if (m != length(n)) 
        stop("objects of different length")
    if (is.null(spChk)) 
        spChk <- get.spChkOption()
    if (spChk && !chkIDs(x, nb)) 
        stop("Check of data and neighbour ID integrity failed")
    if (spChk && !chkIDs(n, nb)) 
        stop("Check of data and neighbour ID integrity failed")
    if (!is.numeric(x)) 
        stop(paste(deparse(substitute(x)), "is not a numeric vector"))
    if (!is.numeric(n)) 
        stop(paste(deparse(substitute(n)), "is not a numeric vector"))
    if (any(is.na(x))) 
        stop("NA in at risk population")
    if (any(is.na(n))) 
        stop("NA in cases")
    if (any(x < 0)) 
        stop("negative risk population")
    if (any(n < 0)) 
        stop("negative number of cases")
    lw <- nb2listw(include.self(nb), style="B", zero.policy=zero.policy)
    r <- n/x
    localcas <- lag.listw(lw, n, zero.policy = zero.policy)
    localpop <- lag.listw(lw, x, zero.policy = zero.policy)
    localmeancas <- lag.listw(nb2listw(include.self(nb), style="W",
        zero.policy=zero.policy), n, zero.policy = zero.policy)

    m <- localcas/localpop
    localC <- lag.listw(lw, (x * (r - m)^2), zero.policy = zero.policy)
    a <- (localC/localpop) - (m/localmeancas)
    a[a < 0] <- 0
    est <- m + (r - m) * (a / (a + (m/x)))

    res <- data.frame(raw=r, est=est)
    attr(res, "parameters") <- list(a=a, m=m)
    res
}
