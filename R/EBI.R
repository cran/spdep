# Copyright 2002-2003 by Roger Bivand 
#

EBImoran <- function(z, listw, nn, S0, zero.policy=FALSE) {
	zm <- mean(z)
	zz <- sum((z - zm)^2)
	lz <- lag.listw(listw, z, zero.policy=zero.policy)
	EBI <- (nn / S0) * ((t(z) %*% lz) / zz)
	res <- EBI
	res
}

EBImoran.mc <- function(n, x, listw, nsim, zero.policy=FALSE,
	alternative="greater", spChk=NULL) {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if(!is.numeric(n)) stop(paste(deparse(substitute(n)),
		"is not a numeric vector"))	
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(x))) stop("NA in at risk population")
	if (any(is.na(n))) stop("NA in cases")
	m <- length(listw$neighbours)
	if (m != length(x)) stop("objects of different length")
	if (m != length(n)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	if (spChk && !chkIDs(n, listw))
		stop("Check of data and weights ID integrity failed")
	if(nsim > gamma(m+1))
		stop("nsim too large for this number of observations")
	S0 <- Szero(listw)
	p <- n/x
	nsum <- sum(n)
	xsum <- sum(x)
	b <- nsum/xsum
	s2 <- sum(x*(((p-b)^2)/xsum))
	a <- s2 - (b/(xsum/m))
	v <- a + (b/x)
	v[v < 0] <- b/x
	z <- (p-b)/sqrt(v)
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) res[i] <- EBImoran(sample(z), listw, m, S0,
	    zero.policy)
	res[nsim+1] <- EBImoran(z, listw, m, S0, zero.policy)
	rankres <- rank(res)
	zrank <- rankres[length(res)]
	diff <- nsim - zrank
	diff <- ifelse(diff > 0, diff, 0)
        pval <- (diff + 1)/(nsim+1)
	if (alternative == "less") pval <- 1 - pval
	else if (alternative == "two.sided") pval <- 2 * pval
	statistic <- res[nsim+1]
	names(statistic) <- "statistic"
	parameter <- zrank
	names(parameter) <- "observed rank"
	method <- "Monte-Carlo simulation of Empirical Bayes Index"
	data.name <- paste("cases: ", deparse(substitute(n)),
	    ", risk population: ", deparse(substitute(x)), "\nweights: ",
	    deparse(substitute(listw)), "\nnumber of simulations + 1: ",
	    nsim+1, "\n", sep="")
	lres <- list(statistic=statistic, parameter=parameter,
	    p.value=pval, alternative=alternative, method=method, 
	    data.name=data.name, res=res, z=z)
	class(lres) <- c("htest", "mc.sim")
	lres
}

