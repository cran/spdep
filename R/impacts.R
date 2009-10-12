# Copyright 2009 by Roger Bivand

trW <- function(W, m=100, p=50, type="mult") {
# returns traces
    n <- dim(W)[1]
    iW <- W
    tr <- numeric(m)
    if (type == "mult") {
        for (i in 1:m) {
            tr[i] <- sum(diag(iW))
            iW <- W %*% iW
        }
    } else if (type == "MC") {
        x <- matrix(rnorm(n*p), nrow=n, ncol=p)
        xx <- x
        for (i in 1:m) {
            xx <- W %*% xx
            tr[i] <- sum(apply(x * as.matrix(xx), 2,  function(y) sum(y)/p))
# mean replaced by sum(y)/p 091012, 0.4-47
        }
        tr[1] <- 0.0
        tr[2] <- sum(t(W) * W)
    } else stop("unknown type")
    tr
}

impacts.sarlm <- function(obj, ..., tr, R=NULL, listw=NULL, useHESS=NULL,
  tol=1e-6, empirical=FALSE) {
    if (obj$type == "error")
        stop("impact measures only for lag and spatial Durbin models")
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
    rho <- obj$rho
    beta <- obj$coefficients
    s2 <- obj$s2
    usingHESS <- NULL
    if (!is.null(R)) {
        resvar <- obj$resvar
        usingHESS <- FALSE
        irho <- 2
        drop2beta <- 1:2
        if (is.logical(resvar)) {
            fdHess <- obj$fdHess
            if (is.logical(fdHess)) 
                stop("coefficient covariance matrix not available")
            usingHESS <- TRUE
            if (!obj$insert) {
                irho <- 1
                drop2beta <- 1
            }
        }
        if (!is.null(useHESS) && useHESS) {
            fdHess <- obj$fdHess
            if (is.logical(fdHess)) 
                stop("Hessian matrix not available")
            usingHESS <- TRUE
            if (!obj$insert) {
                irho <- 1
                drop2beta <- 1
            }
        }
        interval <- obj$interval
        if (is.null(interval)) interval <- c(-1,0.999)
    }
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (obj$type == "lag") {
      if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
      } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
      }
    } else if (obj$type == "mixed") {
      if (iicept) {
        b1 <- beta[-icept]
      } else {
        b1 <- beta
      }
      p <- length(b1)
      if (p %% 2 != 0) stop("non-matched coefficient pairs")
      P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
      bnames <- names(b1[1:(p/2)])
    }
    n <- length(obj$fitted.values)
    if (is.null(listw)) {
        q <- length(tr)-1
        g <- rho^(0:q)
        T <- matrix(c(1, tr[-(q+1)]/n), nrow=1)
        if (obj$type == "mixed") {
            T <- rbind(T, tr/n)
        }
        res <- lagImpacts(T, g, P)
        if (!is.null(R)) {
            if (usingHESS && !obj$insert) {
                mu <- c(rho, beta)
                samples <- mvrnorm(n=R, mu=mu, Sigma=fdHess, tol=tol,
                    empirical=empirical)
            } else {
                mu <- c(s2, rho, beta)
                if (usingHESS) {
                    samples <- mvrnorm(n=R, mu=mu, Sigma=fdHess, tol=tol,
                        empirical=empirical)
                } else {
                    samples <- mvrnorm(n=R, mu=mu, Sigma=resvar, tol=tol,
                        empirical=empirical)
                }
            }
            check <- ((samples[,irho] > interval[1]) & 
                (samples[,irho] < interval[2]))
            if (any(!check)) samples <- samples[check,]
            processSample <- function(x, irho, drop2beta) {
                g <- x[irho]^(0:q)
                beta <- x[-drop2beta]
                if (obj$type == "lag") {
                  if (iicept) {
                    P <- matrix(beta[-icept], ncol=1)
                  } else {
                    P <- matrix(beta, ncol=1)
                  }
                } else if (obj$type == "mixed") {
                    if (iicept) {
                      b1 <- beta[-icept]
                    } else {
                      b1 <- beta
                    }
                    p <- length(b1)
                    if (p %% 2 != 0) stop("non-matched coefficient pairs")
                    P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
                }
                lagImpacts(T, g, P)
            }
            sres <- apply(samples, 1, processSample, irho=irho,
                drop2beta=drop2beta)
            direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
            indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
            total <- as.mcmc(t(sapply(sres, function(x) x$total)))
            colnames(direct) <- bnames
            colnames(indirect) <- bnames
            colnames(total) <- bnames
            res <- list(res=res, sres=list(direct=direct,
                indirect=indirect, total=total))
        }
        attr(res, "method") <- "trace"
    } else {
        SW <- invIrW(listw, rho)
        if (obj$type == "lag") res <- lagImpactsExact(SW, P, n)
        else if (obj$type == "mixed") res <- mixedImpactsExact(SW, P, n, listw)
        if (!is.null(R)) {
            if (usingHESS && !obj$insert) {
                mu <- c(rho, beta)
                samples <- mvrnorm(n=R, mu=mu, Sigma=fdHess, tol=tol,
                    empirical=empirical)
            } else {
                mu <- c(s2, rho, beta)
                if (usingHESS) {
                    samples <- mvrnorm(n=R, mu=mu, Sigma=fdHess, tol=tol,
                        empirical=empirical)
                } else {
                    samples <- mvrnorm(n=R, mu=mu, Sigma=resvar, tol=tol,
                        empirical=empirical)
                }
            }
            check <- ((samples[,irho] > interval[1]) & 
                (samples[,irho] < interval[2]))
            if (any(!check)) samples <- samples[check,]
            processXSample <- function(x, drop2beta) {
                beta <- x[-drop2beta]
                if (obj$type == "lag") {
                    if (iicept) {
                      P <- matrix(beta[-icept], ncol=1)
                    } else {
                      P <- matrix(beta, ncol=1)
                    }
                    return(lagImpactsExact(SW, P, n))
                } else if (obj$type == "mixed") {
                    if (iicept) {
                        b1 <- beta[-icept]
                    } else {
                        b1 <- beta
                    }
                    P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
                    return(mixedImpactsExact(SW, P, n, listw))
                }
            }
            sres <- apply(samples, 1, processXSample, drop2beta=drop2beta)
            direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
            indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
            total <- as.mcmc(t(sapply(sres, function(x) x$total)))
            colnames(direct) <- bnames
            colnames(indirect) <- bnames
            colnames(total) <- bnames
            res <- list(res=res, sres=list(direct=direct,
                indirect=indirect, total=total))
        }
        attr(res, "method") <- "exact"
    }
    attr(res, "useHESS") <- usingHESS
    attr(res, "insert") <- obj$insert
    attr(res, "type") <- obj$type
    attr(res, "bnames") <- bnames
    class(res) <- "sarlmImpact"
    res
}

lagImpacts <- function(T, g, P) {
    PT <- P %*% T
    direct <- apply(apply(PT, 1, function(x) x*g), 2, sum)
    total <- c(apply(P, 1, sum) * sum(g))
    indirect <- total - direct
    names(direct) <- names(total)
    list(direct=direct, indirect=indirect, total=total)
}

lagImpactsExact <- function(SW, P, n) {
    direct <- sapply(P, function(x) sum(diag(x*SW))/n)
    total <- sapply(P, function(x) sum(x*SW)/n)
    indirect <- total - direct
    list(direct=direct, indirect=indirect, total=total)
}

mixedImpactsExact <- function(SW, P, n, listw) {
    p <- dim(P)[1]
    direct <- numeric(p)
    total <- numeric(p)
    W <- listw2mat(listw)
    for (i in 1:p) {
        SWr <- SW %*% (P[i,1]*diag(n) + P[i,2]*W)
        direct[i] <- sum(diag(SWr))/n
        total[i] <- sum(SWr)/n
    }
    indirect <- total - direct
    list(direct=direct, indirect=indirect, total=total)
}

sarlmImpactMat <- function(x) {
    if (is.null(x$res)) {
        direct <- x$direct
        indirect <- x$indirect
        total <- x$total
    } else {
        direct <- x$res$direct
        indirect <- x$res$indirect
        total <- x$res$total
    }
    mat <- cbind(direct, indirect, total)
    colnames(mat) <- c("Direct", "Indirect", "Total")
    rownames(mat) <- attr(x, "bnames")
    mat
}


print.sarlmImpact <- function(x, ...) {
    mat <- sarlmImpactMat(x)
    cat("Impact measures (", attr(x, "type"), ", ", attr(x, "method"), "):\n", sep="")
    
    print(mat)
    invisible(x)
}

summary.sarlmImpact <- function(object, ..., zstats=FALSE, short=FALSE) {
    if (is.null(object$sres)) stop("summary method unavailable")
    direct_sum <- summary(object$sres$direct)
    indirect_sum <- summary(object$sres$indirect)
    total_sum <- summary(object$sres$total)
    lres <- list(direct_sum=direct_sum, indirect_sum=indirect_sum,
        total_sum=total_sum)
    res <- c(object, lres)
    if (zstats) {
        zmat <- sapply(lres, function(x) x$statistics[,1]/x$statistics[,2])
        colnames(zmat) <- c("Direct", "Indirect", "Total")
        pzmat <- 2*(1-pnorm(abs(zmat)))
        res <- c(res, list(zmat=zmat, pzmat=pzmat))
    }
    attr(res, "useHESS") <- attr(object, "useHESS")
    attr(res, "bnames") <- attr(object, "bnames")
    attr(res, "method") <- attr(object, "method")
    attr(res, "insert") <- attr(object, "insert")
    attr(res, "type") <- attr(object, "type")
    attr(res, "short") <- short
    class(res) <- "summary.sarlmImpact"
    res
}

print.summary.sarlmImpact <- function(x, ...) {
    mat <- sarlmImpactMat(x)
    cat("Impact measures (", attr(x, "type"), ", ", attr(x, "method"),
        "):\n", sep="")
    print(mat)
    cat("========================================================\n")
    tp <- ifelse(attr(x, "useHESS"), ifelse(attr(x, "insert"), "mixed Hessian approximation", "numerical Hessian approximation"), "asymptotic")
    cat("Simulation results (", tp, " variance matrix):\n", sep="")
    if (!attr(x, "short")) {
        cat("Direct:\n")
        print(x$direct_sum)
        cat("========================================================\n")
        cat("Indirect:\n")
        print(x$indirect_sum)
        cat("========================================================\n")
        cat("Total:\n")
        print(x$total_sum)
    }
    if (!is.null(x$zmat)) {
        cat("========================================================\n")
        cat("Simulated z-values:\n")
        mat <- x$zmat
        rownames(mat) <- attr(x, "bnames")
        print(mat)
        cat("\nSimulated p-values:\n")
        xx <- apply(x$pzmat, 2, format.pval)
        rownames(xx) <- attr(x, "bnames")
        print(xx, quote=FALSE)
    }
    invisible(x)
}

plot.sarlmImpact <- function(x, ..., choice="direct", trace=FALSE,
    density=TRUE) {
    if (is.null(x$sres)) stop("plot method unavailable")
    plot(x$sres[[choice]], trace=trace, density=density, sub=choice)
    invisible(x)
}

HPDinterval.sarlmImpact <- function(obj, prob = 0.95, ..., choice="direct") {
    if (is.null(obj$sres)) stop("HPDinterval method unavailable")
    res <- HPDinterval(obj$sres[[choice]], prob=prob)
    res
}

