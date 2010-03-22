# Copyright 2009 by Roger Bivand

trW <- function(W, m=100, p=50, type="mult") {
# returns traces
    timings <- list()
    .ptime_start <- proc.time()
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
    timings[["make_traces"]] <- proc.time() - .ptime_start
    attr(tr, "timings") <- do.call("rbind", timings)[, c(1, 3)]
    attr(tr, "type") <- type
    tr
}

impacts <- function(obj, ...)
    UseMethod("impacts", obj)

impacts.stsls <- function(obj, ..., tr=NULL, R=NULL, listw=NULL,
  tol=1e-6, empirical=FALSE, Q=NULL) {
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
    rho <- obj$coefficients[1]
    beta <- obj$coefficients[-1]
    icept <- grep("(Intercept)", names(beta))
    iicept <- length(icept) > 0
    if (iicept) {
        P <- matrix(beta[-icept], ncol=1)
        bnames <- names(beta[-icept])
    } else {
        P <- matrix(beta, ncol=1)
        bnames <- names(beta)
    }
    p <- length(beta)
    n <- length(obj$residuals)
    mu <- c(rho, beta)
    Sigma <- obj$var
    irho <- 1
    drop2beta <- 1
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu, Sigma=Sigma,
        irho=irho, drop2beta=drop2beta, bnames=bnames, interval=NULL,
        type="lag", tr=tr, R=R, listw=listw, tol=tol, empirical=empirical,
        Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "iClass") <- class(obj)
    res
}

intImpacts <- function(rho, beta, P, n, mu, Sigma, irho, drop2beta, bnames,
    interval, type, tr, R, listw, tol, empirical, Q, icept, iicept, p,
    mess=FALSE) {
    if (is.null(listw) && is.null(tr))
        stop("either tr or listw must be given")
    timings <- list()
    .ptime_start <- proc.time()
    if (is.null(listw)) {

        lagImpacts <- function(T, g, P) {
            PT <- P %*% T
            direct <- apply(apply(PT, 1, function(x) x*g), 2, sum)
            total <- c(apply(P, 1, sum) * sum(g))
            indirect <- total - direct
            names(direct) <- names(total)
            list(direct=direct, indirect=indirect, total=total)
        }

        lagDistrImpacts <- function(T, g, P, q=10) {
            PT <- P %*% T
            direct <- apply(PT, 1, function(x) x * g)[1:q, ]
            total <- t(sapply(g, function(x) apply(P, 1, sum)*x))[1:q, ]
            indirect <- total - direct
            list(direct=direct, indirect=indirect, total=total)
        }

        q <- length(tr)-1
        g <- rho^(0:q)
        T <- matrix(c(1, tr[-(q+1)]/n), nrow=1)
        if (type == "mixed") {
            T <- rbind(T, tr/n)
        }
        res <- lagImpacts(T, g, P)
        if (!is.null(Q)) {
            if (!is.numeric(Q) || length(Q) > 1) stop("Invalid Q argument")
            if (Q > length(tr)) stop("Q larger than length of tr")
            Qres <- lagDistrImpacts(T, g, P, q=as.integer(Q))
            attr(res, "Qres") <- Qres
        }
        timings[["trace_impacts"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        if (!is.null(R)) {
            samples <- mvrnorm(n=R, mu=mu, Sigma=Sigma, tol=tol,
                empirical=empirical)
            if (mess) samples[,irho] <- 1 - exp(samples[,irho])
            if (!is.null(interval)) {
                check <- ((samples[,irho] > interval[1]) & 
                    (samples[,irho] < interval[2]))
                if (any(!check)) samples <- samples[check,]
            }
            timings[["impacts_samples"]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
            processSample <- function(x, irho, drop2beta) {
                g <- x[irho]^(0:q)
                beta <- x[-drop2beta]
                if (type == "lag") {
                  if (iicept) {
                    P <- matrix(beta[-icept], ncol=1)
                  } else {
                    P <- matrix(beta, ncol=1)
                  }
                } else if (type == "mixed") {
                    if (iicept) {
                      b1 <- beta[-icept]
                    } else {
                      b1 <- beta
                    }
                    p <- length(b1)
                    if (p %% 2 != 0) stop("non-matched coefficient pairs")
                    P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
                }
                res <- lagImpacts(T, g, P)
                if (!is.null(Q)) {
                    Qres <- lagDistrImpacts(T, g, P, q=as.integer(Q))
                    attr(res, "Qres") <- Qres
                }
                res
            }
            CL <- get("cl", env = .spdepOptions)
            if (!is.null(CL) && length(CL) > 1) {
                require(snow)
                l_sp <- lapply(splitIndices(nrow(samples), length(CL)), 
		    function(i) samples[i,])
		clusterExport_l <- function(CL, list) {
                    gets <- function(n, v) {
                        assign(n, v, env = .GlobalEnv)
                        NULL
                    }
                    for (name in list) {
                        clusterCall(CL, gets, name, get(name))
                    }
		}
		clusterExport_l(CL, list("processSample", "irho", "drop2beta",
                    "Q", "T", "lagImpacts", "lagDistrImpacts", "icept",
                    "iicept", "type"))
                timings[["cluster_setup"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                lsres <- parLapply(CL, l_sp, function(sp) apply(sp, 1, 
                    processSample, irho=irho, drop2beta=drop2beta))
		clusterEvalQ(CL, rm(list=c("processSample", "irho", "drop2beta",
                    "Q", "T", "lagImpacts", "lagDistrImpacts", "icept",
                    "iicept", "type")))
                sres <- do.call("c", lsres)

            } else {
                sres <- apply(samples, 1, processSample, irho=irho,
                    drop2beta=drop2beta)
            }
            timings[["process_samples"]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
            direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
            indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
            total <- as.mcmc(t(sapply(sres, function(x) x$total)))
            colnames(direct) <- bnames
            colnames(indirect) <- bnames
            colnames(total) <- bnames
            ssres <- list(direct=direct, indirect=indirect, total=total)
            if (!is.null(Q)) {
                Qdirect <- as.mcmc(t(sapply(sres, function(x)
                    attr(x, "Qres")$direct)))
                Qindirect <- as.mcmc(t(sapply(sres, function(x) 
                    attr(x, "Qres")$indirect)))
                Qtotal <- as.mcmc(t(sapply(sres, function(x) 
                    attr(x, "Qres")$total)))
                Qnames <- c(sapply(bnames, function(x) 
                    paste(x, 1:Q, sep="__Q")))
                colnames(Qdirect) <- Qnames
                colnames(Qindirect) <- Qnames
                colnames(Qtotal) <- Qnames
                Qmcmc <- list(direct=Qdirect, indirect=Qindirect, total=Qtotal)
                attr(ssres, "Qmcmc") <- Qmcmc
            }
            timings[["postprocess_samples"]] <- proc.time() - .ptime_start
            res <- list(res=res, sres=ssres)
        }
        attr(res, "method") <- "trace"
    } else {
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
        SW <- invIrW(listw, rho)
        if (type == "lag") res <- lagImpactsExact(SW, P, n)
        else if (type == "mixed") res <- mixedImpactsExact(SW, P, n, listw)
        timings[["weights_impacts"]] <- proc.time() - .ptime_start
        .ptime_start <- proc.time()
        if (!is.null(R)) {
            samples <- mvrnorm(n=R, mu=mu, Sigma=Sigma, tol=tol,
                empirical=empirical)
            if (!is.null(interval)) {
                check <- ((samples[,irho] > interval[1]) & 
                    (samples[,irho] < interval[2]))
                if (any(!check)) samples <- samples[check,]
            }
            timings[["impacts_samples"]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
            processXSample <- function(x, drop2beta) {
                beta <- x[-drop2beta]
                if (type == "lag") {
                    if (iicept) {
                      P <- matrix(beta[-icept], ncol=1)
                    } else {
                      P <- matrix(beta, ncol=1)
                    }
                    return(lagImpactsExact(SW, P, n))
                } else if (type == "mixed") {
                    if (iicept) {
                        b1 <- beta[-icept]
                    } else {
                        b1 <- beta
                    }
                    P <- cbind(b1[1:(p/2)], b1[((p/2)+1):p])
                    return(mixedImpactsExact(SW, P, n, listw))
                }
            }
            CL <- get("cl", env = .spdepOptions)
            if (!is.null(CL) && length(CL) > 1) {
                require(snow)
                l_sp <- lapply(splitIndices(nrow(samples), length(CL)), 
		    function(i) samples[i,])
		clusterExport_l <- function(CL, list) {
                    gets <- function(n, v) {
                        assign(n, v, env = .GlobalEnv)
                        NULL
                    }
                    for (name in list) {
                        clusterCall(CL, gets, name, get(name))
                    }
		}

		clusterExport_l(CL, list("processXSample", "irho", "drop2beta",
                    "SW", "lagImpactsExact", "mixedImpactsExact", "icept",
                    "iicept", "type", "listw"))

                timings[["cluster_setup"]] <- proc.time() - .ptime_start
                .ptime_start <- proc.time()
                lsres <- parLapply(CL, l_sp, function(sp) apply(sp, 1, 
                    processXSample, drop2beta=drop2beta))
		clusterEvalQ(CL, rm(list=c("processXSample", "drop2beta",
                    "SW", "lagImpactsExact", "mixedImpactsExact", "icept",
                    "iicept", "type", "listw")))
                sres <- do.call("c", lsres)

            } else {
                sres <- apply(samples, 1, processXSample,
                    drop2beta=drop2beta)
            }
            timings[["process_samples"]] <- proc.time() - .ptime_start
            .ptime_start <- proc.time()
            direct <- as.mcmc(t(sapply(sres, function(x) x$direct)))
            indirect <- as.mcmc(t(sapply(sres, function(x) x$indirect)))
            total <- as.mcmc(t(sapply(sres, function(x) x$total)))
            colnames(direct) <- bnames
            colnames(indirect) <- bnames
            colnames(total) <- bnames
            timings[["postprocess_samples"]] <- proc.time() - .ptime_start
            res <- list(res=res, sres=list(direct=direct,
                indirect=indirect, total=total))
        }
        attr(res, "method") <- "exact"
    }
    attr(res, "type") <- type
    attr(res, "bnames") <- bnames
    attr(res, "haveQ") <- !is.null(Q)
    attr(res, "timings") <- do.call("rbind", timings)[, c(1,3)]
    class(res) <- "lagImpact"
    res
}

impacts.lagmess <- function(obj, ..., tr=NULL, R=NULL, listw=NULL, 
  tol=1e-6, empirical=FALSE, Q=NULL) {
    stopifnot(!is.null(obj$mixedHess))
    stop("method not yet available")
}

impacts.sarlm <- function(obj, ..., tr=NULL, R=NULL, listw=NULL, useHESS=NULL,
  tol=1e-6, empirical=FALSE, Q=NULL) {
    if (obj$type == "error")
        stop("impact measures only for lag and spatial Durbin models")
    if (!is.null(obj$listw_style) && obj$listw_style != "W") 
        stop("Only row-standardised weights supported")
    rho <- obj$rho
    beta <- obj$coefficients
    s2 <- obj$s2
    usingHESS <- NULL
    iNsert <- obj$insert
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
            if (!iNsert) {
                irho <- 1
                drop2beta <- 1
            }
        }
        if (!is.null(useHESS) && useHESS) {
            fdHess <- obj$fdHess
            if (is.logical(fdHess)) 
                stop("Hessian matrix not available")
            usingHESS <- TRUE
            if (!iNsert) {
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
      p <- length(beta)
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
    n <- length(obj$residuals)
    mu <- NULL
    Sigma <- NULL
    if (!is.null(R)) {
        if (usingHESS && !iNsert) {
            mu <- c(rho, beta)
            Sigma <- fdHess
        } else {
            mu <- c(s2, rho, beta)
            if (usingHESS) {
                Sigma <- fdHess
            } else {
                Sigma <- resvar
            }
        }
    }
    res <- intImpacts(rho=rho, beta=beta, P=P, n=n, mu=mu, Sigma=Sigma,
        irho=irho, drop2beta=drop2beta, bnames=bnames, interval=interval,
        type=obj$type, tr=tr, R=R, listw=listw, tol=tol, empirical=empirical,
        Q=Q, icept=icept, iicept=iicept, p=p)
    attr(res, "useHESS") <- usingHESS
    attr(res, "insert") <- iNsert
    attr(res, "iClass") <- class(obj)
    res
}

lagImpactMat <- function(x, reportQ=NULL) {
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
    if (!is.null(reportQ) && reportQ) {
        if (is.null(x$res)) {
            Qobj <- attr(x, "Qres")
        } else {
            Qobj <- attr(x$res, "Qres")
        }
        if (is.null(Qobj)) warning("No impact components to report")
        else {
            colnames(Qobj$direct) <- attr(x, "bnames")
            colnames(Qobj$indirect) <- attr(x, "bnames")
            colnames(Qobj$total) <- attr(x, "bnames")
            attr(mat, "Qobj") <- Qobj
        }
    }
    mat
}


print.lagImpact <- function(x, ..., reportQ=NULL) {
    mat <- lagImpactMat(x, reportQ=reportQ)
    Qobj <- attr(mat, "Qobj")
    cat("Impact measures (", attr(x, "type"), ", ", attr(x, "method"), "):\n", sep="")
    attr(mat, "Qobj") <- NULL
    print(mat)
    if (!is.null(reportQ) && reportQ) {
        if (is.null(Qobj)) warning("No impact components to report")
        else {
            cat("=================================\nImpact components\n")
            print(Qobj)
        }
    }
    invisible(x)
}

summary.lagImpact <- function(object, ..., zstats=FALSE, short=FALSE, reportQ=NULL) {
    if (is.null(object$sres)) stop("summary method unavailable")
    direct_sum <- summary(object$sres$direct)
    indirect_sum <- summary(object$sres$indirect)
    total_sum <- summary(object$sres$total)
    Qmcmc <- NULL
    if (!is.null(attr(object$sres, "Qmcmc")) && !is.null(reportQ) && reportQ) {
        Qdirect_sum <- summary(attr(object$sres, "Qmcmc")$direct)
        Qindirect_sum <- summary(attr(object$sres, "Qmcmc")$indirect)
        Qtotal_sum <- summary(attr(object$sres, "Qmcmc")$total)
        Qmcmc <- list(Qdirect_sum=Qdirect_sum, Qindirect_sum=Qindirect_sum,
           Qtotal_sum=Qtotal_sum)
    }
    lres <- list(direct_sum=direct_sum, indirect_sum=indirect_sum,
        total_sum=total_sum)
    res <- c(object, lres, Qmcmc)
    if (zstats) {
        zmat <- sapply(lres, function(x) x$statistics[,1]/x$statistics[,2])
        colnames(zmat) <- c("Direct", "Indirect", "Total")
        pzmat <- 2*(1-pnorm(abs(zmat)))
        res <- c(res, list(zmat=zmat, pzmat=pzmat))
        if (!is.null(Qmcmc) && !is.null(reportQ) && reportQ) {
            Qzmats <- lapply(Qmcmc, function(x) {
                Qm <- matrix(x$statistics[,1]/x$statistics[,2],
                    ncol=length(attr(object, "bnames")))
                colnames(Qm) <- attr(object, "bnames")
                Qm
            })
            names(Qzmats) <- c("Direct", "Indirect", "Total")
            Qpzmats <- lapply(Qzmats, function(x) 2*(1-pnorm(abs(x))))
            res <- c(res, list(Qzmats=Qzmats, Qpzmats=Qpzmats))
        }
    }
    attr(res, "useHESS") <- attr(object, "useHESS")
    attr(res, "bnames") <- attr(object, "bnames")
    attr(res, "method") <- attr(object, "method")
    attr(res, "insert") <- attr(object, "insert")
    attr(res, "type") <- attr(object, "type")
    attr(res, "short") <- short
    attr(res, "reportQ") <- reportQ
    tp <- NULL
    if (attr(object, "iClass") == "sarlm") tp <- ifelse(attr(object,
       "useHESS"), ifelse(attr(object, "insert"),
       "mixed Hessian approximation", "numerical Hessian approximation"),
       "asymptotic")
    else if (attr(object, "iClass") == "stsls") tp <- "asymptotic IV"
    attr(res, "tp") <- tp
    class(res) <- "summary.lagImpact"
    res
}

print.summary.lagImpact <- function(x, ...) {
    reportQ <- attr(x, "reportQ")
    mat <- lagImpactMat(x, reportQ)
    Qobj <- attr(mat, "Qobj")
    attr(mat, "Qobj") <- NULL
    cat("Impact measures (", attr(x, "type"), ", ", attr(x, "method"),
        "):\n", sep="")
    print(mat)
    if (!is.null(reportQ) && reportQ) {
        if (is.null(Qobj)) warning("No impact components to report")
        else {
            cat("=================================\nImpact components\n")
            print(Qobj)
        }
    }
    cat("========================================================\n")

    cat("Simulation results (", attr(x, "tp"), " variance matrix):\n", sep="")
    if (!attr(x, "short")) {
        cat("Direct:\n")
        print(x$direct_sum)
        cat("========================================================\n")
        cat("Indirect:\n")
        print(x$indirect_sum)
        cat("========================================================\n")
        cat("Total:\n")
        print(x$total_sum)
        if (!is.null(reportQ) && reportQ && !is.null(x$Qdirect_sum)) {
            cat("========================================================\n")
            cat("Direct impact components:\n")
            print(x$Qdirect_sum)
            cat("========================================================\n")
            cat("Indirect impact components:\n")
            print(x$Qindirect_sum)
            cat("========================================================\n")
            cat("Total impact components:\n")
            print(x$Qtotal_sum)
        }
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
        if (!is.null(x$Qzmats)) {
            cat("========================================================\n")
            cat("Simulated impact components z-values:\n")
            print(x$Qzmats)
            cat("\nSimulated impact components p-values:\n")
            xx <- lapply(x$Qpzmats, function(y) apply(y, 2, format.pval))
            print(xx, quote=FALSE)
        }
    }
    invisible(x)
}

plot.lagImpact <- function(x, ..., choice="direct", trace=FALSE,
    density=TRUE) {
    if (is.null(x$sres)) stop("plot method unavailable")
    plot(x$sres[[choice]], trace=trace, density=density, sub=choice)
    invisible(x)
}

HPDinterval.lagImpact <- function(obj, prob = 0.95, ..., choice="direct") {
    if (is.null(obj$sres)) stop("HPDinterval method unavailable")
    res <- HPDinterval(obj$sres[[choice]], prob=prob)
    res
}

