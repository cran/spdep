# Copyright 2001-2 by Roger Bivand 
#

lm.morantest <- function(model, listw, zero.policy=FALSE, 
	    alternative = "greater") {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(class(model) != "lm") stop(paste(deparse(substitute(model)),
		"not an lm object"))
 	N <- length(listw$neighbours)
	u <- as.vector(residuals(model))
	if (N != length(u)) 
            stop("objects of different length")
	listw.U <- listw2U(listw)

	S0 <- sum(unlist(listw.U$weights))
	S1 <- 0.5 * sum((2*unlist(listw.U$weights))^2)
	lu <- lag.listw(listw.U, u, zero.policy=zero.policy)
	I <- (N/S0) * ((t(u) %*% lu) / (t(u) %*% u))
	p <- model$rank
	p1 <- 1:p
	XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
	X <- model.matrix(terms(model), model.frame(model))
# Cliff/Ord 1981, p. 203
	Z <- lag.listw(listw.U, X, zero.policy=zero.policy)
	C1 <- t(X) %*% Z
	trA <- -(sum(diag(XtXinv %*% C1)))
	EI <- ((N * trA) / ((N-p) * S0))
	C2 <- t(Z) %*% Z
	C3 <- XtXinv %*% C1
	trA2 <- sum(diag(C3 %*% C3))
	trB <- sum(diag(4*(XtXinv %*% C2)))
	VI <- (((N*N)/((S0*S0)*(N-p)*(N-p+2))) *
		(S1 + 2*trA2 - trB - ((2*(trA^2))/(N-p))))
	ZI <- (I - EI) / sqrt(VI)
    	if (alternative == "two.sided") pv <- 2 * (1 - pnorm(ZI))
    	else if (alternative == "greater")
	        pv <- pnorm(ZI, lower.tail=FALSE)
    	else pv <- pnorm(ZI)
    	statistic <- ZI
    	attr(statistic, "names") <- "Saddlepoint approximation"
    	p.value <- pv
    	estimate <- c(I, EI, VI)
    	attr(estimate, "names") <- c("Observed Moran's I", "Expectation",
	    "Variance")
    	method <- "Global Moran's I for regression residuals"
    	data.name <- paste("model:", deparse(model$call),
    	    "\nweights:", deparse(substitute(listw)), "\n")
    	res <- list(statistic = statistic, p.value = p.value,
	       estimate = estimate, method = method,
		alternative = alternative, data.name = data.name)
	class(res) <- "htest"
    	res
}
