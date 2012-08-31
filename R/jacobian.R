# Copyright 2010-2011 by Roger Bivand

# Chebyshev approximation setup and run functions
cheb_setup <- function(env, q=5, which=1) {
    if (which == 1) {
        W <- as(as_dgRMatrix_listw(get("listw", envir=env)), "CsparseMatrix")
    } else {
        W <- as(as_dgRMatrix_listw(get("listw2", envir=env)), "CsparseMatrix")
    }
# W a CSparseMatrix object
# q order
    n <- nrow(W)
    IM <- .symDiagonal(n)
    T <- vector(mode="list", length=(q+1))
    T[[1]] <- IM
    T[[2]] <- W
    trT <- numeric(q+1)
    trT[1] <- n
    trT[2] <- 0
    if (q > 1) {
        for (k in 3:(q+1)) {
            T[[k]] <- 2*(W %*% T[[(k-1)]]) - T[[(k-2)]]
            trT[k] <- sum(diag(T[[k]]))
        }
    }
    if (which == 1) {
      assign("trT", trT, envir=env)
      assign("W", W, envir=env)
    } else {
      assign("trT2", trT, envir=env)
      assign("W2", W, envir=env)
    }
    assign("method", "Chebyshev", envir=env)
    invisible(NULL)
}

cheb_ldet <- function(alpha, env, which=1) {
# trT output from cheb_setup()
# alpha spatial coefficient
    if (which == 1) {
        trT <- get("trT", envir=env)
    } else {
        trT <- get("trT2", envir=env)
    }
    cheb_in <- function(alpha, j, q) {
        res <- (2/(q+1))
        x <- 0.0
        for (k in 1:(q+1)) {
            x <- x + log(((1 - (alpha*cos((pi*(k - 0.5))/(q + 1)))))) * 
                cos((pi*(j - 1)*(k - 0.5))/(q + 1))
        }
        res <- res * x
        res
    }
    q <- length(trT)-1L
    n <- trT[1]
    C1 <- cheb_in(alpha, j=1, q)
    x <- 0.0
    for (j in 1:(q+1)) {
        x <- x + (cheb_in(alpha, j=j, q)*trT[j])
    }
    x <- x - (n/2)*C1
    x
}

# MC approximation setup and run functions
mcdet_setup <- function(env, p=16, m=30, which=1) {
        if (which == 1) {
          W <- as(as_dgRMatrix_listw(get("listw", envir=env)), "CsparseMatrix")
        } else {
          W <- as(as_dgRMatrix_listw(get("listw2", envir=env)), "CsparseMatrix")
        }
# W a CSparseMatrix object
# p, m given in papers
        n <- dim(W)[1]
	x <- matrix(rnorm(n*p), nrow=n, ncol=p)
        int1  <- vector(mode="list", length=m)
	xx <- x
# 111110 set first two traces
        td2 <- sum(diag(W %*% W))
        for (k in 1:m) {
            xx <- W %*% xx
            if (k == 1) int1[[k]] <- rep(0, p)
            else if (k == 2) int1[[k]] <- rep(td2, p)
            else int1[[k]] <- apply(x * as.matrix(xx), 2, sum)
        }
        int2 <- apply(x * x, 2, sum)
        clx <- list(m=m, p=p, n=n, int1=int1, int2=int2)
        if (which == 1) {
            assign("clx", clx, envir=env)
            assign("W", W, envir=env)
        } else {
            assign("clx2", clx, envir=env)
            assign("W2", W, envir=env)
        }
        assign("method", "MC", envir=env)
        invisible(NULL)
}

mcdet_ldet <- function(alpha, env, which=1) {
# clx output from mcdet_setup()
# alpha spatial coefficient
        if (which == 1) {
            clx <- get("clx", envir=env)
        } else {
            clx <- get("clx2", envir=env)
        }
	vk <- numeric(length=clx$p)
	for (k in 1:clx$m) {
		vk <- clx$n*(alpha^k)*(clx$int1[[k]]/k) + vk
	}
	v <- -c(as.matrix(vk/clx$int2))
	res <- mean(v)
        attr(res, "sd") <- sd(v)/sqrt(clx$p)
        res
}

eigen_setup <- function(env, which=1) {
    if (get("verbose", envir=env))
       cat("Computing eigenvalues ...\n")
    if (which == 1) {
	if (get("listw", envir=env)$style %in% c("W", "S") && 
            get("can.sim", envir=env)) {
            eig <- eigen(similar.listw_Matrix(get("listw", envir=env)),
                only.values=TRUE)$value
	    assign("similar", TRUE, envir=env)
	} else eig <- eigenw(get("listw", envir=env))
# modified 110414 RSB
	if (is.complex(eig)) eig.range <- 1/range(Re(eig[which(Im(eig) == 0)]))
	else eig.range <- 1/range(eig)
        assign("eig", eig, envir=env)
        assign("eig.range", eig.range, envir=env)
    } else {
	if (get("listw2", envir=env)$style %in% c("W", "S") && 
            get("can.sim2", envir=env)) {
            eig <- eigen(similar.listw_Matrix(get("listw2", envir=env)),
                only.values=TRUE)$value
	    assign("similar2", TRUE, envir=env)
	} else eig <- eigenw(get("listw2", envir=env))
        assign("eig2", eig, envir=env)
    }
    if (get("verbose", envir=env)) cat("\n")
    assign("method", "eigen", envir=env)
    invisible(NULL)
}

do_ldet <- function(coef, env, which=1) {
    method <- get("method", envir=env)
    if (get("family", envir=env) == "SMA") {
        ldet <- eigen_sma_ldet(coef, env, which=which)
    } else {
        switch(method,
           eigen = {ldet <- eigen_ldet(coef, env, which=which)},
           spam = {ldet <- spam_ldet(coef, env, which=which)},
           spam_update = {ldet <- spam_update_ldet(coef, env, which=which)},
           Matrix = {ldet <- Matrix_ldet(coef, env, which=which)},
           Matrix_J = {ldet <- Matrix_J_ldet(coef, env, which=which)},
           Chebyshev = {ldet <- cheb_ldet(coef, env, which=which)},
           MC = {ldet <- mcdet_ldet(coef, env, which=which)},
           LU = {ldet <- LU_ldet(coef, env, which=which)},
           moments = {ldet <- moments_ldet(coef, env, which=which)},
           SE_classic = {ldet <- SE_classic_ldet(coef, env, which=which)},
           SE_whichMin = {ldet <- SE_whichMin_ldet(coef, env, which=which)},
           SE_interp = {ldet <- SE_interp_ldet(coef, env, which=which)},
           stop("...\n\nUnknown method\n"))
    }
    ldet
}

eigen_sma_ldet <- function(coef, env, which=1) {
    eig <- get("eig", envir=env)
# modified 110414 RSB
    if (is.complex(eig)) det <- Re(sum(log(1/(1 + coef * eig))))
    else det <- sum(log(1/(1 + coef * eig)))
    det
}

eigen_ldet <- function(coef, env, which=1) {
    if (which == 1) {
        eig <- get("eig", envir=env)
    } else {
        eig <- get("eig2", envir=env)
    }
# modified 110414 RSB
    if (is.complex(eig)) 
        det <- Re(sum(log(1 - coef * eig)))
    else det <- sum(log(1 - coef * eig))
    det
}

spam_setup <- function(env, pivot="MMD", which=1) {
    if (!require(spam)) stop("spam not available")
    if (which == 1) {
        if (get("listw", envir=env)$style %in% c("W", "S") &&
            get("can.sim", envir=env)) {
	      csrw <- listw2U_spam(similar.listw_spam(get("listw", envir=env)))
	      assign("similar", TRUE, envir=env)
	} else csrw <- as.spam.listw(get("listw", envir=env))
        assign("csrw", csrw, envir=env)
    } else {
        if (get("listw2", envir=env)$style %in% c("W", "S") &&
            get("can.sim2", envir=env)) {
	      csrw <- listw2U_spam(similar.listw_spam(get("listw2", envir=env)))
	      assign("similar2", TRUE, envir=env)
	} else csrw <- as.spam.listw(get("listw2", envir=env))
        assign("csrw2", csrw, envir=env)
    }
    n <- get("n", envir=env)
    I <- diag.spam(1, n, n)
    assign("I", I, envir=env)
    assign("pivot", pivot, envir=env)
    assign("method", "spam", envir=env)
    invisible(NULL)
}

spam_ldet <- function(coef, env, which=1) {
    if (!require(spam)) stop("spam not available")
    if (which == 1) {
        csrw <- get("csrw", envir=env)
    } else {
        csrw <- get("csrw2", envir=env)
    }
    I <- get("I", envir=env)
    pivot <- get("pivot", envir=env)
    J1 <- try(determinant(chol((I - coef * csrw), pivot=pivot),
        logarithm=TRUE)$modulus, silent=TRUE)
    if (class(J1) == "try-error") {
        Jacobian <- NA
    } else {
        Jacobian <- 2*J1
    }
    Jacobian
}

spam_update_setup <- function(env, in_coef=0.1, pivot="MMD", which=1) {
    if (!require(spam)) stop("spam not available")
    if (which == 1) {
        if (get("listw", envir=env)$style %in% c("W", "S") &&
            get("can.sim", envir=env)) {
	      csrw <- listw2U_spam(similar.listw_spam(get("listw", envir=env)))
	      assign("similar", TRUE, envir=env)
	} else csrw <- as.spam.listw(get("listw", envir=env))
        assign("csrw", csrw, envir=env)
    } else {
        if (get("listw2", envir=env)$style %in% c("W", "S") &&
            get("can.sim2", envir=env)) {
	      csrw <- listw2U_spam(similar.listw_spam(get("listw2", envir=env)))
	      assign("similar2", TRUE, envir=env)
	} else csrw <- as.spam.listw(get("listw2", envir=env))
        assign("csrw2", csrw, envir=env)
    }
    n <- get("n", envir=env)
    I <- diag.spam(1, n, n)
    assign("I", I, envir=env)
    csrwchol <- chol((I - in_coef * csrw), pivot=pivot)
    if (which == 1) {
        assign("csrwchol", csrwchol, envir=env)
    } else {
        assign("csrwchol2", csrwchol, envir=env)
    }
    assign("method", "spam_update", envir=env)
    invisible(NULL)
}

spam_update_ldet <- function(coef, env, which=1) {
    if (!require(spam)) stop("spam not available")
    if (which == 1) {
        csrw <- get("csrw", envir=env)
        cchol <- get("csrwchol", envir=env)
    } else {
        csrw <- get("csrw2", envir=env)
        cchol <- get("csrwchol2", envir=env)
    }
    I <- get("I", envir=env)
    if (abs(coef) < .Machine$double.eps^(0.5)) {
        Jacobian <- 0.0
    } else {
        J1 <- try(determinant(update(cchol, (I - coef * csrw)),
            logarithm=TRUE)$modulus, silent=TRUE)
        if (class(J1) == "try-error") {
            Jacobian <- NA
        } else {
            Jacobian <- 2*J1
        }
    }
    Jacobian
}

Matrix_setup <- function(env, Imult, super=as.logical(NA), which=1) {
    if (which == 1) {
        if (get("listw", envir=env)$style %in% c("W", "S") && 
            get("can.sim", envir=env)) {
	    csrw <- listw2U_Matrix(similar.listw_Matrix(get("listw", 
                envir=env)))
	    assign("similar", TRUE, envir=env)
	} else csrw <- as_dsTMatrix_listw(get("listw", envir=env))
	csrw <- as(csrw, "CsparseMatrix")
        nW <- - csrw
	pChol <- Cholesky(csrw, super=super, Imult = Imult)
	nChol <- Cholesky(nW, super=super, Imult = Imult)
        assign("csrw", csrw, envir=env)
        assign("nW", nW, envir=env)
        assign("pChol", pChol, envir=env)
        assign("nChol", nChol, envir=env)
    } else {
        if (get("listw2", envir=env)$style %in% c("W", "S") && 
            get("can.sim2", envir=env)) {
	    csrw <- listw2U_Matrix(similar.listw_Matrix(get("listw2", 
                envir=env)))
	    assign("similar2", TRUE, envir=env)
	} else csrw <- as_dsTMatrix_listw(get("listw2", envir=env))
	csrw <- as(csrw, "CsparseMatrix")
        nW <- - csrw
	pChol <- Cholesky(csrw, super=super, Imult = Imult)
	nChol <- Cholesky(nW, super=super, Imult = Imult)
        assign("csrw2", csrw, envir=env)
        assign("nW2", nW, envir=env)
        assign("pChol2", pChol, envir=env)
        assign("nChol2", nChol, envir=env)
    }
    .f <- if (package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1
    assign(".f", .f, envir=env)
    assign("method", "Matrix", envir=env)
    invisible(NULL)
}

Matrix_ldet <- function(coef, env, which=1) {
    if (which == 1) {
        csrw <- get("csrw", envir=env)
        nW <- get("nW", envir=env)
        pChol <- get("pChol", envir=env)
        nChol <- get("nChol", envir=env)
    } else {
        csrw <- get("csrw2", envir=env)
        nW <- get("nW2", envir=env)
        pChol <- get("pChol2", envir=env)
        nChol <- get("nChol2", envir=env)
    }
    a <- -.Machine$double.eps^(1/2)
    b <- .Machine$double.eps^(1/2)
    n <- get("n", envir=env)
    .f <- get(".f", envir=env)

    Jacobian <- ifelse(coef > b, n * log(coef) +
            (.f * c(determinant(update(nChol, nW, 1/coef))$modulus)),
            ifelse(coef < a, n* log(-(coef)) + 
            (.f * c(determinant(update(pChol, csrw, 1/(-coef)))$modulus)),
            0.0))
    Jacobian
}

LU_setup <- function(env, which=1) {
    if (which == 1) {
        W <- as(as_dgRMatrix_listw(get("listw", envir=env)), "CsparseMatrix")
        assign("W", W, envir=env)
    } else {
        W <- as(as_dgRMatrix_listw(get("listw2", envir=env)), "CsparseMatrix")
        assign("W2", W, envir=env)
    }
    I <- as_dsCMatrix_I(get("n", envir=env))
    assign("I", I, envir=env)
    assign("method", "LU", envir=env)
    invisible(NULL)
}

LU_ldet <- function(coef, env, which=1) {
    I <- get("I", envir=env)
    if (which == 1) {
        W <- get("W", envir=env)
    } else {
        W <- get("W2", envir=env)
    }
    LU <- lu(I - coef * W)
    dU <- abs(diag(slot(LU, "U")))
    ldet <- sum(log(dU))
    ldet
}

Matrix_J_setup <- function(env, super=FALSE, which=1) {
    if (which == 1) {
        if (get("listw", envir=env)$style %in% c("W", "S") && 
            get("can.sim", envir=env)) {
            csrw <- listw2U_Matrix(similar.listw_Matrix(get("listw", 
                envir=env)))
	    assign("similar", TRUE, envir=env)
        } else csrw <- as_dsTMatrix_listw(get("listw", envir=env))
        csrw <- as(csrw, "CsparseMatrix")
        assign("csrw", csrw, envir=env)
    } else {
        if (get("listw2", envir=env)$style %in% c("W", "S") && 
            get("can.sim2", envir=env)) {
	    csrw <- listw2U_Matrix(similar.listw_Matrix(get("listw2", 
                envir=env)))
	    assign("similar2", TRUE, envir=env)
	} else csrw <- as_dsTMatrix_listw(get("listw2", envir=env))
	csrw <- as(csrw, "CsparseMatrix")
        assign("csrw2", csrw, envir=env)
    }
    I <- as_dsCMatrix_I(get("n", envir=env))
    assign("I", I, envir=env)
    .f <- if (package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1
    assign(".f", .f, envir=env)
    assign("super", super, envir=env)
    assign("method", "Matrix_J", envir=env)
    invisible(NULL)
}

Matrix_J_ldet <- function(coef, env, which=1) {
    I <- get("I", envir=env)
    super <- get("super", envir=env)
    if (which == 1) {
        csrw <- get("csrw", envir=env)
    } else {
        csrw <- get("csrw2", envir=env)
    }
    .f <- get(".f", envir=env)
    cch <- Cholesky((I - coef * csrw), super=super)
    Jacobian <- .f * determinant(cch, logarithm = TRUE)$modulus
    Jacobian
}


Rmrho <- function(Omega, m, rho, n, trunc=FALSE) {
    Om <- Omega[m]
    Om1 <- Omega[m-1]
    Om_e <- Omega[m]/Omega[m-2]
    Om_o <- Omega[m-1]/Omega[m-3]
    res <- 0
    rhoj <- rho^m
    Om_ej <- Om_e^m
    Om_oj <- Om_o^m
    for (j in m:n) {
        if ((j %% 2) == 0) {
            inc <- ((1/j)*rhoj)*Om*(Om_ej)
        } else { 
            inc <- ((1/j)*rhoj)*Om1*(Om_oj)
        }
        if (!is.finite(inc)) break
        if (abs(inc) < .Machine$double.eps && trunc) break
        res <- res + inc
        rhoj <- rhoj*rho
        Om_ej <- Om_ej*Om_e
        Om_oj <- Om_oj*Om_o
    }
    attr(res, "j") <- j
    res
}

ldetMoments <- function(Omega, rho, n, correct=TRUE, trunc=FALSE) {
    m <- length(Omega)
    Rm <- 0
    attr(Rm, "j") <- as.integer(NA)
    if (correct) Rm <- Rmrho(Omega, m, rho, n, trunc)
    res <- 0
    rhoj <- rho
    for (j in seq(along=Omega)) {
        res <- res + (1/j)*rhoj*Omega[j]
        rhoj <- rhoj*rho
    }
    res <- - res - Rm
    attr(res, "j") <- attr(Rm, "j")
    res
}

moments_setup <- function(env, trs=NULL, m, p, type="MC", correct=TRUE,
    trunc=TRUE, which=1) {
    if (which == 1) {
        if (is.null(trs)) {
            if (get("listw", envir=env)$style %in% c("W", "S") && 
                get("can.sim", envir=env)) {
                csrw <- listw2U_Matrix(similar.listw_Matrix(get("listw",
                    envir=env)))
	        assign("similar", TRUE, envir=env)
            } else csrw <- as_dgRMatrix_listw(get("listw", envir=env))
            csrw <- as(csrw, "CsparseMatrix")
            trs <- trW(csrw, m=m, p=p, type=type)
        }
        assign("trs1", trs, envir=env)
    } else {
        if (is.null(trs)) {
            if (get("listw2", envir=env)$style %in% c("W", "S") && 
                get("can.sim2", envir=env)) {
                csrw <- listw2U_Matrix(similar.listw_Matrix(get("listw2",
                    envir=env)))
	        assign("similar2", TRUE, envir=env)
            } else csrw <- as_dgRMatrix_listw(get("listw2", envir=env))
            csrw <- as(csrw, "CsparseMatrix")
            trs <- trW(csrw, m=m, p=p, type=type)
        }
        assign("trs2", trs, envir=env)
    }
    assign("correct", correct, envir=env)
    assign("trunc", trunc, envir=env)
    assign("method", "moments", envir=env)
    invisible(NULL)
}

moments_ldet <- function(x, env, which=1) {
    n <- get("n", envir=env)
    if (which == 1) {
        trs <- get("trs1", envir=env)
    } else {
        trs <- get("trs2", envir=env)
    }
    correct <- get("correct", envir=env)
    trunc <- get("trunc", envir=env)
    Jacobian <- ldetMoments(trs, x, n, correct, trunc)
    Jacobian
}

SE_classic_setup <- function(env, SE_method="LU", p=16, m=30, nrho=200,
  interpn=2000, interval=c(-1,0.999), SElndet=NULL, which=1) {
  stopifnot(require(splines))

  if (is.null(SElndet)) {
    SE_setup_intern(env, SE_method=SE_method, p=p, m=m, nrho=nrho,
      interval=interval, which=which)
    assign("SE_method", SE_method, envir=env)
    if (which == 1) {
      detval <- get("detval1", envir=env)
    } else if (which == 2) {
      detval <- get("detval2", envir=env)
    }
    fit <- interpSpline(detval[,1], detval[,2])
    rho <- seq(interval[1], interval[2], length.out=interpn)
    detval <- matrix(unlist(predict(fit, rho)), ncol=2)
  } else {
    stopifnot(is.matrix(SElndet))
    stopifnot(ncol(SElndet) == 2)
    detval <- SElndet
    assign("SE_method", "precomputed", envir=env)
  }
  assign("method", "SE_classic", envir=env)
  if (which == 1) {
    assign("detval1", detval, envir=env)
  } else if (which == 2) {
    assign("detval2", detval, envir=env)
  }
  assign("intern_classic", data.frame(), envir=env)
  invisible(NULL)
 
}

SE_setup_intern <- function(env, SE_method="LU", p=16, m=30, nrho=100,
  interval=c(-1,0.999), which=1) {

  switch(SE_method,
    LU = {tull <- LU_setup(env, which=which)},
    MC = {tull <- mcdet_setup(env, p=16, m=30, which=which)}, 
    stop("...\n\nUnknown SE_method\n"))

    rho <- seq(interval[1], interval[2], length.out=nrho)

    ldets <- sapply(rho, function(r) do_ldet(r, env, which=which))
    detval <- cbind(rho, ldets)
    
    if (which == 1) {
        assign("detval1", detval, envir=env)
    } else if (which == 2) {
        assign("detval2", detval, envir=env)
    }

    invisible(NULL)
}

SE_classic_ldet <- function(x, env, which=1) {
    if (which == 1) {
        detval <- get("detval1", envir=env)
    } else if (which == 2) {
        detval <- get("detval2", envir=env)
    }
    res <- SE_classic(x, detval)
    intern_attr <- attr(res, "intern")
    intern_attr$rho0 <- x
    intern_attr$rho1 <- res[1]
    intern_attr <- as.data.frame(intern_attr)
    iC <- rbind(get("intern_classic", envir=env), intern_attr)
    assign("intern_classic", iC, envir=env)
    res[2]
}


SE_classic <- function(rho, detval) {

  gsize = detval[2, 1] - detval[1, 1]
  i1 = which(detval[, 1] <= rho + gsize)
  i2 = which(detval[, 1] <= rho - gsize)
  i1 = max(i1)
  i2 = max(i2)
  index0 <- (i1+i2)/2
  index = round(index0)
#cat("index", index, "i1", i1, "i2", i2, "ind", ((i1+i2)/2), "\n")
  if (index < 1 || index > dim(detval)[1]) stop("index out of bounds")

  res <- detval[index, ]
  attr(res, "intern") <- list(i1=i1, i2=i2, index0=index0, index=index)
  res

}

SE_whichMin_setup <- function(env, SE_method="LU", p=16, m=30, nrho=200,
  interpn=2000, interval=c(-1,0.999), SElndet=NULL, which=1) {
  stopifnot(require(splines))

  if (is.null(SElndet)) {
    SE_setup_intern(env, SE_method=SE_method, p=p, m=m, nrho=nrho,
      interval=interval, which=which)
    assign("SE_method", SE_method, envir=env)
    if (which == 1) {
      detval <- get("detval1", envir=env)
    } else if (which == 2) {
      detval <- get("detval2", envir=env)
    }
  
    fit <- interpSpline(detval[,1], detval[,2])
    rho <- seq(interval[1], interval[2], length.out=interpn)
    detval <- matrix(unlist(predict(fit, rho)), ncol=2)
  } else {
    stopifnot(is.matrix(SElndet))
    stopifnot(ncol(SElndet) == 2)
    detval <- SElndet
    assign("SE_method", "precomputed", envir=env)
  }
  assign("method", "SE_whichMin", envir=env)
  if (which == 1) {
    assign("detval1", detval, envir=env)
  } else if (which == 2) {
    assign("detval2", detval, envir=env)
  }

  invisible(NULL)
}

SE_whichMin_ldet <- function(x, env, which=1) {
    if (which == 1) {
        detval <- get("detval1", envir=env)
    } else if (which == 2) {
        detval <- get("detval2", envir=env)
    }
    SE_whichMin(x, detval)[2]
}


SE_whichMin <- function(rho, detval) {

  gsize = detval[2, 1] - detval[1, 1]
  i1 = which(detval[, 1] <= rho + gsize)
  i2 = which(detval[, 1] <= rho - gsize)
  i1 = max(i1)
  i2 = max(i2)
  i12 <- i1:i2
  mi12 <- which.min((detval[i12, 1]-rho)^2)
  index <- i12[mi12]
#cat("index", index, "i1", i1, "i2", i2, "\n")
  if (index < 1 || index > dim(detval)[1]) stop("index out of bounds")

  detval[index, ]

}

SE_interp_setup <- function(env, SE_method="LU", p=16, m=30, nrho=200,
  interval=c(-1,0.999), which=1) {

  stopifnot(require(splines))

  SE_setup_intern(env, SE_method=SE_method, p=p, m=m, nrho=nrho,
    interval=interval, which=which)
  assign("method", "SE_interp", envir=env)
  assign("SE_method", SE_method, envir=env)
  if (which == 1) {
    detval <- get("detval1", envir=env)
  } else if (which == 2) {
    detval <- get("detval2", envir=env)
  }
  fit <- interpSpline(detval[,1], detval[,2])
  if (which == 1) {
    assign("fit1", fit, envir=env)
  } else if (which == 2) {
    assign("fit2", fit, envir=env)
  }

  invisible(NULL)
}

SE_interp_ldet <- function(x, env, which=1) {
    if (which == 1) {
        fit <- get("fit1", envir=env)
    } else if (which == 2) {
        fit <- get("fit2", envir=env)
    }
    SE_interp(x, fit)[2]
}


# setup generate detval using method # for # rho, fit spline model

SE_interp <- function(rho, fit) {

  res <- predict(fit, rho)
  unname(unlist(res))
}

