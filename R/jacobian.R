# Copyright 2010 by Roger Bivand

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
    q <- length(trT)-1
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
        for (k in 1:m) {
            xx <- W %*% xx
            int1[[k]] <- apply(x * as.matrix(xx), 2, sum)
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
	v <- c(as.matrix(vk/clx$int2))
	-mean(v)
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
	if (is.complex(eig)) eig.range <- 1/range(Re(eig))
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
           stop("...\n\nUnknown method\n"))
    }
    ldet
}

eigen_sma_ldet <- function(coef, env, which=1) {
    eig <- get("eig", envir=env)
    if (is.complex(eig)) det <- sum(log(Re(1/(1 + coef * eig))))
    else det <- sum(log(1/(1 + coef * eig)))
    det
}

eigen_ldet <- function(coef, env, which=1) {
    if (which == 1) {
        eig <- get("eig", envir=env)
    } else {
        eig <- get("eig2", envir=env)
    }
    if (is.complex(eig)) 
        det <- sum(log(1 - coef * Re(eig)))
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

Matrix_setup <- function(env, Imult, super, which=1) {
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

Matrix_J_setup <- function(env, which=1) {
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
    assign("method", "Matrix_J", envir=env)
    invisible(NULL)
}

Matrix_J_ldet <- function(coef, env, which=1) {
    I <- get("I", envir=env)
    if (which == 1) {
        csrw <- get("csrw", envir=env)
    } else {
        csrw <- get("csrw2", envir=env)
    }
    Jacobian <- determinant(I - coef * csrw, logarithm = TRUE)$modulus
    Jacobian
}

