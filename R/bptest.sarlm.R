# Copyright 2004 by Roger Bivand (original taken from bptest() in the lmtest
# package, Copyright (C) 2001 Torsten Hothorn and Achim Zeileis and released
# under GNU General Public License, Version 2.
#

bptest.sarlm <- function (object, studentize = TRUE) 
{
    if(!inherits(object, "sarlm")) stop("not sarlm object")
    Z <- model.matrix(object$lm.target)
    k <- ncol(Z)
    n <- nrow(Z)
    resi <- residuals(object$lm.target)
    sigma2 <- sum(resi^2)/n
    if (studentize) {
        w <- resi^2 - sigma2
        fv <- lm.fit(Z, w)$fitted
        bp <- n * sum(fv^2)/sum(w^2)
        method <- "studentized Breusch-Pagan test"
    }
    else {
        f <- resi^2/sigma2 - 1
        fv <- lm.fit(Z, f)$fitted
        bp <- 0.5 * sum(fv^2)
        method <- "Breusch-Pagan test"
    }
    names(bp) <- "BP"
    df <- k - 1
    names(df) <- "df"
    RVAL <- list(statistic = bp, parameter = df, method = method, 
        p.value = 1 - pchisq(bp, df))
    class(RVAL) <- "htest"
    return(RVAL)
}

