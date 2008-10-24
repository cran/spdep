aple.plot <- function(x, listw, do.plot=TRUE, ...) {
    pre <- preAple(x=x, listw=listw)
    W2e <- eigen(pre$W2)
    SQRTW2 <- W2e$vectors %*% (diag(W2e$values^(0.5)) %*% t(W2e$vectors))
    X <- drop(SQRTW2 %*% x)
    NSQRTW2 <- W2e$vectors %*% (diag(W2e$values^(-0.5)) %*% t(W2e$vectors))
    Y <- drop(NSQRTW2 %*% pre$WU %*% x)
    if (do.plot) {
        plot(X, Y, ...)
    }
    list(X=X, Y=Y)
}

localAple <- function(x, listw) {
    aplepl <- aple.plot(x, listw, do.plot=FALSE)
    res <- (length(aplepl$X) * aplepl$Y * aplepl$X) / crossprod(aplepl$X)
    res
}
