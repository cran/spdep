# Copyright 2001-2 by Roger Bivand 
#

moran <- function(x, listw, n, S0, zero.policy=FALSE) {
	z <- scale(x, scale=FALSE)
	zz <- sum(z^2)
	K <- (n*sum(z^4))/(zz^2)
	lz <- lag.listw(listw, z, zero.policy)
	I <- (n / S0) * ((t(z) %*% lz) / zz)
	res <- list(I=I, K=K)
	res
}

moran.test <- function(x, listw, randomisation=TRUE, zero.policy=FALSE,
	alternative="greater") {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	wc <- spweights.constants(listw)
	S02 <- wc$S0*wc$S0
	res <- moran(x, listw, wc$n, wc$S0, zero.policy)
	I <- res$I
	K <- res$K
	EI <- (-1) / wc$n1
	if(randomisation) {
		VI <- n*(wc$S1*(wc$nn - 3*n + 3) - n*wc$S2 + 3*S02)
		tmp <- K*(wc$S1*(wc$nn - n) - 2*n*wc$S2 + 6*S02)
		VI <- (VI - tmp) / (wc$n1*wc$n2*wc$n3*S02)
		VI <- VI - EI^2
	} else {
		VI <- (wc$nn*wc$S1 - n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
		VI <- VI - EI^2
	}
	ZI <- (I - EI) / sqrt(VI)
	statistic <- ZI
	names(statistic) <- "Moran I statistic standard deviate"
        if (alternative == "two.sided") PrI <- 2 * pnorm(ZI)
        else if (alternative == "greater")
            PrI <- pnorm(ZI, lower.tail=FALSE)
        else PrI <- pnorm(ZI)
	vec <- c(I, EI, VI)
	names(vec) <- c("Moran I statistic", "Expectation", "Variance")
	method <- paste("Moran's I test under", ifelse(randomisation,
	    "randomisation", "normality"))
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\n")
	res <- list(statistic=statistic, p.value=PrI, estimate=vec, 
	    alternative=alternative, method=method, data.name=data.name)
	class(res) <- "htest"
	res
}

moran.mc <- function(x, listw, nsim, zero.policy=FALSE) {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	S0 <- Szero(listw)
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) res[i] <- moran(sample(x), listw, n, S0,
	    zero.policy)$I
	res[nsim+1] <- moran(x, listw, n, S0, zero.policy)$I
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	statistic <- res[nsim+1]
	names(statistic) <- "Moran's I statistic"
	parameter <- xrank
	names(parameter) <- "rank of observed I"
	method <- "Monte-Carlo simulation of Moran's I"
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\nnumber of simulations + 1:",
	    nsim+1, "\n")
	lres <- list(statistic=statistic, parameter=parameter,
	    method=method, data.name=data.name, res=res)
	class(lres) <- "htest"
	lres
}


