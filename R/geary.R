# Copyright 2001-2 by Roger Bivand 
#


geary <- function(x, listw, n, n1, S0, zero.policy=FALSE) {
	z <- scale(x, scale=FALSE)
	zz <- sum(z^2)
	K <- (n*sum(z^4))/(zz^2)
	res <- geary.intern(x, listw, n, zero.policy, type="geary")
	C <- (n1 / (2*S0)) * (sum(res) / zz)
	res <- list(C=C, K=K)
	res
}

#geary.intern <- function(x, listw, n, zero.policy, type="geary") {
#	res <- as.numeric(rep(0, n))
#	cardnb <- card(listw$neighbours)
#	if (type == "geary") f <- function(diff) (diff^2)
#	else if (type == "sokal") f <- function(diff) (abs(diff))
#	else stop("type unknown")
#	for (i in 1:n) {
#		if (cardnb[i] == 0) {
#			if (zero.policy) res[i] <- 0
#			else res[i] <- NA
#		} else {
#			res[i] <- sum(listw$weights[[i]] * 
#				f(x[i] - x[listw$neighbours[[i]]]))
#		}
#	}
#	res
#}

geary.intern <- function(x, listw, n, zero.policy, type="geary") {
	cardnb <- card(listw$neighbours)
	if (type == "geary") ft <- TRUE
	else if (type == "sokal") ft <- FALSE
	else stop("type unknown")
	res <- .Call("gearyw", listw$neighbours, listw$weights,
		as.numeric(x), as.integer(cardnb),
		as.logical(zero.policy), as.logical(ft), PACKAGE="spdep")
	if (any(is.na(res))) warning("NAs in lagged values")
	invisible(res)
}

geary.test <- function(x, listw, randomisation=TRUE, zero.policy=FALSE,
    alternative="less") {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	wc <- spweights.constants(listw, zero.policy)
	S02 <- wc$S0*wc$S0
	res <- geary(x, listw, wc$n, wc$n1, wc$S0, zero.policy)
	C <- res$C
	if (is.na(C)) stop("NAs generated in geary - check zero.policy")
	K <- res$K
	EC <- 1
	if(randomisation) {
		VC <- (wc$n1*wc$S1*(wc$nn - 3*n + 3 - K*wc$n1))
		VC <- VC - ((1/4) * (wc$n1*wc$S2*(wc$nn + 3*n - 6 - 
			K*(wc$nn - n + 2))))
		VC <- VC + (S02*(wc$nn - 3 - K*(wc$n1^2)))
		VC <- VC / (n*wc$n2*wc$n3*S02)
	} else {
		VC <- ((2*wc$S1 + wc$S2)*wc$n1 - 4*S02) / (2*(n + 1)*S02)
	}
	ZC <- (C - EC) / sqrt(VC)
	statistic <- ZC
	names(statistic) <- "Geary C statistic standard deviate"
        if (alternative == "two.sided") PrC <- 2 * pnorm(ZC)
        else if (alternative == "greater")
            PrC <- pnorm(ZC, lower.tail=FALSE)
        else PrC <- pnorm(ZC)
	vec <- c(C, EC, VC)
	names(vec) <- c("Geary C statistic", "Expectation", "Variance")
	method <- paste("Geary's C test under", ifelse(randomisation,
	    "randomisation", "normality"))
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\n")
	res <- list(statistic=statistic, p.value=PrC, estimate=vec, 
	    alternative=alternative, method=method, data.name=data.name)
	class(res) <- "htest"
	res
}

geary.mc <- function(x, listw, nsim, zero.policy=FALSE,
	alternative="less") {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if(nsim > gamma(n+1)) stop("nsim too large for this n")
	wc <- spweights.constants(listw, zero.policy)
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) res[i] <- geary(sample(x), listw, n, wc$n1, wc$S0,
	    zero.policy)$C
	res[nsim+1] <- geary(x, listw, n, wc$n1, wc$S0, zero.policy)$C
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	diff <- nsim - xrank
	diff <- ifelse(diff > 0, diff, 0)
        pval <- (diff + 1)/(nsim+1)
	if (alternative == "less") pval <- 1 - pval
	else if (alternative == "two.sided") pval <- 2 * pval
	statistic <- res[nsim+1]
	names(statistic) <- "statistic"
	parameter <- xrank
	names(parameter) <- "observed rank"
	method <- "Monte-Carlo simulation of Geary's C"
	data.name <- paste(deparse(substitute(x)), "\nweights:",
	    deparse(substitute(listw)), "\nnumber of simulations + 1:",
	    nsim+1, "\n")
	lres <- list(statistic=statistic, parameter=parameter,
	    p.value=pval, alternative=alternative, method=method, 
	    data.name=data.name, res=res)
	class(lres) <- c("htest", "mc.sim")
	lres
}

