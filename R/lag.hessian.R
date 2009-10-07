f_laglm_eig <- function(coefs, y, X, yW, n, eig) {
    rho <- coefs[1]
    beta <- coefs[-1]
    res <- (y - rho * yW) - X %*% beta
    SSE <- sum(res^2)
    s2 <- SSE/n
    if (is.complex(eig)) 
        det <- sum(log(1 - rho * Re(eig)))
    else det <- sum(log(1 - rho * eig))
    ret <- (det - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
   ret
}

getVmat_eig <- function(coefs, y, X, yW, n, eig, s2, trs, tol.solve=1.0e-10,
    optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_eig, y=y, X=X, yW=yW, n=n, eig=eig,
            method="BFGS", hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_laglm_eig, y, X, yW, n, eig)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asy(coefs, y, X, yW, n, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

f_laglm_Matrix <- function(coefs, y, X, yW, n, W, I, nW, nChol, pChol) {
    rho <- coefs[1]
    beta <- coefs[-1]
    res <- (y - rho * yW) - X %*% beta
    SSE <- sum(res^2)
    s2 <- SSE/n
    a <- -.Machine$double.eps^(1/2)
    b <- .Machine$double.eps^(1/2)
    .f <- if (package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1

    Jacobian <- ifelse(rho > b, n * log(rho) +
        (.f * c(determinant(update(nChol, nW, 1/rho))$modulus)),
        ifelse(rho < a, n* log(-(rho)) + (.f * c(determinant(update(pChol,
        W, 1/(-rho)))$modulus)), 0.0))
    ret <- (Jacobian - ((n/2) * log(2 * pi)) - (n/2) * log(s2) - 
        (1/(2 * s2)) * SSE)
   ret
}

getVmat_Matrix <- function(coefs, y, X, yW, n, W, I, nW, nChol, pChol, s2,
    trs, tol.solve=1.0e-10, optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_Matrix, y=y, X=X, yW=yW, n=n,
            W=W, I=I, nW=nW, nChol=nChol, pChol=pChol, method="BFGS",
            hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_laglm_Matrix, y, X, yW, n, W, I, nW, nChol,
            pChol)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asy(coefs, y, X, yW, n, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

f_laglm_spam <- function(coefs, y, X, yW, n, W, I) {
    rho <- coefs[1]
    beta <- coefs[-1]
    res <- (y - rho * yW) - X %*% beta
    SSE <- sum(res^2)
    s2 <- SSE/n
    J1 <- try(determinant((I - rho * W), logarithm=TRUE)$modulus,
        silent=TRUE)
   if (class(J1) == "try-error") {
      	Jacobian <- NA
   } else {
      	Jacobian <- J1
   }
   ret <- (Jacobian - ((n/2)*log(2*pi)) - (n/2)*log(s2) - (1/(2*s2))*SSE)
   ret
}


getVmat_spam <- function(coefs, y, X, yW, n, W, I, s2, trs, tol.solve=1.0e-10,
    optim=FALSE) {
    if (optim) {
        opt <- optim(par=coefs, fn=f_laglm_spam, y=y, X=X, yW=yW, n=n,
            W=W, I=I, method="BFGS", hessian=TRUE)
        mat <- opt$hessian
    } else {
        fd <- fdHess(coefs, f_laglm_spam, y, X, yW, n, W, I)
        mat <- fd$Hessian
    }
    if (!is.null(trs)) {
         mat <- insert_asy(coefs, y, X, yW, n, s2, mat, trs)
    }
    res <- solve(-(mat), tol.solve=tol.solve)
    res
}

trB <- function(rho, tr)  sum(sapply(0:(length(tr)-1),
    function(i) rho^i * tr[i+1]))

insert_asy <- function(coefs, y, X, yW, n, s2, mat, trs) {
    p <- length(coefs)-1
    p2 <- p+2
    omat <- matrix(0, nrow=p2, ncol=p2)
    omat[3:p2, 3:p2] <- -crossprod(X)/s2
    omat[2, 2] <- mat[1, 1]
    omat[2, 3:p2] <- omat[3:p2, 2] <- -c(crossprod(yW, X)/s2)
    omat[1, 1] <- -n/(2*(s2^2))
    omat[1, 2] <- omat[2, 1] <- -trB(coefs[1], trs)/s2
    omat
}

