preAple <- function(x, listw, override_similarity_check=FALSE) {
    stopifnot(isTRUE(all.equal(mean(x), 0.0)))
    if (listw$style %in% c("W", "S") && !override_similarity_check) {
        can.sim <- can.be.simmed(listw)
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

aple <- function(x, listw, override_similarity_check=FALSE) {
    pre <- preAple(x=x, listw=listw,
        override_similarity_check=override_similarity_check)
    res <- inAple(x=x, pre=pre)
    res
}


