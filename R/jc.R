# Copyright 2001-2 by Roger Bivand 
#

joincount <- function(dums, listw) {
	nc <- ncol(dums)
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	res <- as.numeric(rep(0, nc))
	for (lev in 1:nc) {
		for (i in 1:n) {
			xi <- dums[i, lev]
			if (cardnb[i] > 0)
				res[lev] <- res[lev] + (dums[i, lev] *
				sum(dums[listw$neighbours[[i]], lev] *
				listw$weights[[i]]))
		}
	}
	res
}

joincount.test <- function(fx, listw,
	alternative="greater") {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.factor(fx)) stop(paste(deparse(substitute(x)),
		"is not a factor"))
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	wc <- spweights.constants(listw)
	S02 <- wc$S0*wc$S0

	dums <- lm(codes(fx) ~ fx - 1, x=TRUE)$x
	BB <- joincount(dums, listw)
	nBB <- length(BB)
	res <- vector(mode="list", length=nBB)
	tab <- table(fx)
	BB5 <- 0.5 * BB
	ntab <- as.vector(tab)
	Ejc <- (wc$S0*(ntab*(ntab-1))) / (2*n*wc$n1)
	Vjc <- (wc$S1*(ntab*(ntab-1))) / (n*wc$n1)
	Vjc <- Vjc + (((wc$S2 - 2*wc$S1)*ntab*(ntab-1)*(ntab-2)) /
		(n*wc$n1*wc$n2))
	Vjc <- Vjc + (((S02 + wc$S1 - wc$S2)*ntab*(ntab-1)*(ntab-2)*
		(ntab-3)) / (n*wc$n1*wc$n2*wc$n3))
	Vjc <- (0.25 * Vjc) - Ejc^2
	for (i in 1:nBB) {
		estimate <- c(BB5[i], Ejc[i], Vjc[i])
		names(estimate) <- c("Same colour statistic",
			"Expectation", "Variance")
		statistic <- (BB5[i] - Ejc[i]) / sqrt(Vjc[i])
		names(statistic) <- paste("Std. deviate for", names(tab)[i])
		if (alternative == "two.sided") p.value <- 2 * pnorm(statistic)
		else if (alternative == "greater")
			p.value <- pnorm(statistic, lower.tail=FALSE)
		else p.value <- pnorm(statistic)
		method <- "Join count test under nonfree sampling"
		data.name <- paste(deparse(substitute(fx)), "\nweights:",
			deparse(substitute(listw)), "\n")
		res[[i]] <- list(statistic=statistic, p.value=p.value,
			estimate=estimate, method=method,
			alternative=alternative, data.name=data.name)
		class(res[[i]]) <- "htest"
	}
	class(res) <- "jclist"
	res
}

print.jclist <- function(x, ...) {
	for (i in 1:length(x)) print(x[[i]], ...)
	invisible(x)
}

joincount.mc <- function(fx, listw, nsim) {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.factor(fx)) stop(paste(deparse(substitute(fx)),
		"is not a factor"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	dums <- lm(codes(fx) ~ fx - 1, x=TRUE)$x
	nc <- ncol(dums)
	res <- matrix(0, nrow=nsim+1, ncol=nc)
	res[nsim+1,] <- 0.5 * joincount(dums, listw)
	tab <- table(fx)
	for (i in 1:nsim) {
		fxi <- sample(fx)
		dums <- lm(codes(fxi) ~ fxi - 1, x=TRUE)$x
		res[i,] <- 0.5 * joincount(dums, listw)
	}
	rankres <- apply(res, 2, rank)
	xrank <- rankres[nrow(rankres),]
	lres <- vector(mode="list", length=nc)
	for (i in 1:nc) {
		statistic <- res[nrow(res), i]
		names(statistic) <- paste("Join-count statistic for",
			names(tab)[i])
		parameter <- xrank[i]
		names(parameter) <- "rank of observed statistic"
		method <- "Monte-Carlo simulation of join-count statistic"
		data.name <- paste(deparse(substitute(fx)), "\nweights:",
			deparse(substitute(listw)),
			"\nnumber of simulations + 1:", nsim+1, "\n")
		estimate <- c(mean(res[-(nrow(res)), i]),
			var(res[-(nrow(res)), i]))
		names(estimate) <- c("mean of simulation",
			"variance of simulation")
		lres[[i]] <- list(statistic=statistic, parameter=parameter,
			method=method, data.name=data.name,
			estimate=estimate, res=res[,i])
		class(lres[[i]]) <- "htest"
		
	}
	class(lres) <- "jclist"
	lres
}


