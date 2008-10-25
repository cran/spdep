preAple <- function(x, listw) {
    stopifnot(isTRUE(all.equal(mean(x), 0.0)))
    if (listw$style %in% c("W", "S")) {
        can.sim <- spdep:::can.be.simmed(listw)
        eig <- eigenw(similar.listw(listw))
    } else {
        can.sim <- FALSE
        eig <- eigenw(listw)
    }
    if (is.complex(eig)) 
        eig <- Re(eig)
    n <- length(eig)
    corterm <- (crossprod(eig)/n) * Diagonal(n)
    corterm <- as(corterm, "CsparseMatrix")
    W <- as_dgRMatrix_listw(listw)
    W <- as(W, "CsparseMatrix")
    WU <- ((W + t(W))/2)
    W2 <- crossprod(W) + corterm
    res <- list(W=W, corterm=corterm, W2=W2, WU=WU, n=n)
    res
}

inAple <- function(x, pre) {
    xwx <- crossprod(x, (pre$WU %*% x))
    xwwx <- crossprod(x, (pre$W2 %*% x))
    res <- c(as.matrix(xwx/xwwx))
    res
}

aple <- function(x, listw) {
    pre <- preAple(x=x, listw=listw)
    res <- inAple(x=x, pre=pre)
    res
}


