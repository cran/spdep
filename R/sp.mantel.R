# Copyright 2002 by Roger Bivand 
#

sp.mantel.mc <- function(var, listw, nsim, type="moran", zero.policy=FALSE,
	alternative="greater", spChk=NULL) {
	alternative <- match.arg(alternative, c("greater", "less"))
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(var)) stop(paste(deparse(substitute(var)),
		"is not a numeric vector"))
	if(missing(nsim)) stop("nsim must be given")
	n <- length(listw$neighbours)
	if(nsim > gamma(n+1)) stop("nsim too large for this n")
	if (any(is.na(var))) stop("NA in var")
	if (n != length(var)) stop("objects of different length")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(var, listw))
		stop("Check of data and weights ID integrity failed")
	
	listw.U <- listw2U(listw)
	mantel.moran <- function(x, listwU, zero.policy) {
		res <- x * lag.listw(listw.U, x, zero.policy=zero.policy)
		res <- sum(res)
		res
	}
	mantel.geary <- function(x, listwU, zero.policy) {
		res <- geary.intern(x, listwU, length(x), 
			zero.policy=zero.policy, type="geary")
		res <- sum(res)
		res
	} 
	mantel.sokal <- function(x, listwU, zero.policy) {
		res <- geary.intern(x, listwU, length(x), 
			zero.policy=zero.policy, type="sokal")
		res <- sum(res)
		res
	}
	if (type == "moran") f <- mantel.moran
	else if (type == "geary") f <- mantel.geary
	else if (type == "sokal") f <- mantel.sokal
	else stop("unknown type")
	xs <- scale(var)
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) {
		y <- sample(xs)
		res[i] <- f(y, listw.U, zero.policy=zero.policy)
	}
	res[nsim+1] <- f(xs, listw.U, zero.policy=zero.policy)
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	diff <- nsim - xrank
	diff <- ifelse(diff > 0, diff, 0)
        if (alternative == "less") 
        	pval <- punif((diff + 1)/(nsim + 1), lower.tail=FALSE)
    	else if (alternative == "greater") 
        	pval <- punif((diff + 1)/(nsim + 1))
	if (pval < 0 || pval > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	statistic <- res[nsim+1]
	names(statistic) <- "statistic"
	parameter <- xrank
	names(parameter) <- "observed rank"
	method <- paste("Mantel permutation test for", type, "measure")
	data.name <- paste(deparse(substitute(var)), "\nweights:",
	    deparse(substitute(listw)), "\nnumber of simulations + 1:",
	    nsim+1, "\n")
	est <- c(mean(res[1:nsim]), sd(res[1:nsim]))
	names(est) <- c("mean of permutations", "sd of permutations")
	lres <- list(statistic=statistic, parameter=parameter,
	    alternative=alternative, method=method, data.name=data.name, 
	    p.value=pval, res=res, estimate=est)
	class(lres) <- c("htest", "mc.sim")
	lres
}

plot.mc.sim <- function(x, ...) {
	res <- x$res
	xlim <- range(res)
	n <- length(res)
	obs <- res[n]
	res <- res[-n]
	plot(density(res), xlim=xlim, xlab=strsplit(x$data.name, "\n")[[1]][1], 
		main="Density plot of permutation outcomes", sub=x$method)
	abline(v=obs)
}
