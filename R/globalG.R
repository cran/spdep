# Copyright 2002 by Hisaji ONO and Roger Bivand 
#
# General G Statistics
#
#
globalG.test <- function(x, listw, zero.policy=FALSE,
	alternative="greater") {
	if (class(listw) != "listw")
	stop(paste(deparse(substitute(listw)), "is not a listw object"))
	if (listw$style != "B") stop("Only binary weights allowed")
	if (!is.numeric(x))
	stop(paste(deparse(substitute(x)), "is not a numeric vector"))
	if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
	if (any(x < 0.0)) 
		stop(paste("Negative value in ", deparse(substitute(x))))
	n <- length(listw$neighbours)
	if (n != length(x))stop("Different numbers of observations")

	wc <- spweights.constants(listw, zero.policy=zero.policy)
	n1 <- n - 1
	n2 <- n - 2
	n3 <- n - 3
	nn <- n * n
	S0 <- wc$S0
	S1 <- wc$S1
	S2 <- wc$S2
	S02 <- S0*S0
	G <- (t(x) %*% lag.listw(listw, x, zero.policy=zero.policy)) /
		(sum(x %x% x) - (t(x) %*% x))

	E.G <- S0 / (n * n1)

	B0 <- ((nn - 3*n + 3)*S1) - (n*S2) + (3*S02)
	B1 <- -(((nn - n)*S1) - (2*n*S2) + (6*S02))
	B2 <- -((2*n*S1) - ((n+3)*S2) + (6*S02))
	B3 <- (4*n1*S1) - (2*(n+1)*S2) + (8*S02)
	B4 <- S1 - S2 + S02
	sx <- sum(x)
	sx2 <- sum(x^2)
	sx3 <- sum(x^3)
	sx4 <- sum(x^4)

	var.G <- ((B0*(sx2^2) + B1*sx4 + B2*(sx^2)*sx2 + B3*sx*sx3 +
		 B4*(sx^4)) / ((((sx^2) - sx2)^2)*n*n1*n2*n3)) - (E.G^2)

	statistic <- (G - E.G) / sqrt(var.G)
	names(statistic) <- "Getis-Ord global G statistic"
	if (alternative == "two.sided") PrG <- 2 * pnorm(statistic, 
		lower.tail=FALSE)
        else if (alternative == "greater")
            PrG <- pnorm(statistic, lower.tail=FALSE)
        else PrI <- pnorm(statistic)
	vec <- c(G, E.G, var.G)
	names(vec) <- c("Global G statistic", "Expectation", "Variance")
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\n")
	res <- list(statistic=statistic, p.value=PrG, estimate=vec, 
	    alternative=alternative, data.name=data.name)
	class(res) <- "htest"
	res
}

