# Copyright 2001-3 by Roger Bivand 
#

joincount <- function(dums, listw) {
	nc <- which(colSums(dums) > 1)
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	res <- as.numeric(rep(0, ncol(dums)))
	for (lev in nc) {
		res[lev] <- .Call("jcintern", listw$neighbours,
			listw$weights, as.integer(dums[,lev]),
			as.integer(cardnb), PACKAGE="spdep")
	}
	res
}

joincount.test <- function(fx, listw, zero.policy=FALSE,
	alternative="greater", spChk=NULL) {
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.factor(fx)) stop(paste(deparse(substitute(x)),
		"is not a factor"))
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	cards <- card(listw$neighbours)
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(fx, listw))
		stop("Check of data and weights ID integrity failed")
	if (!(alternative %in% c("greater", "less", "two.sided")))
		stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
	wc <- spweights.constants(listw, zero.policy=zero.policy)
	S02 <- wc$S0*wc$S0

	ff <- ~ fx - 1
	dums <- model.matrix(ff, model.frame(ff))
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
		if (alternative == "two.sided") 
			p.value <- 2 * pnorm(-abs(statistic))
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

joincount.mc <- function(fx, listw, nsim, zero.policy=FALSE,
	alternative="greater", spChk=NULL) {
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.factor(fx)) stop(paste(deparse(substitute(fx)),
		"is not a factor"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	cards <- card(listw$neighbours)
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(fx, listw))
		stop("Check of data and weights ID integrity failed")
	if(nsim > gamma(n+1)) stop("nsim too large for this n")
	if (!(alternative %in% c("greater", "less", "two.sided")))
		stop("alternative must be one of: \"greater\", \"less\", or \"two.sided\"")
	ff <- ~ fx - 1
	dums <- model.matrix(ff, model.frame(ff))
	nc <- ncol(dums)
	res <- matrix(0, nrow=nsim+1, ncol=nc)
	res[nsim+1,] <- 0.5 * joincount(dums, listw)
	tab <- table(fx)
	for (i in 1:nsim) {
		fxi <- sample(fx)
		ff <- ~ fxi - 1
		dums <- model.matrix(ff, model.frame(ff))
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
		diff <- nsim - xrank[i]
		diff <- ifelse(diff > 0, diff, 0)
        	pval <- (diff + 1)/(nsim+1)
		if (alternative == "less") pval <- 1 - pval
		else if (alternative == "two.sided") pval <- 2 * pval
		method <- "Monte-Carlo simulation of join-count statistic"
		data.name <- paste(deparse(substitute(fx)), "\nweights:",
			deparse(substitute(listw)),
			"\nnumber of simulations + 1:", nsim+1, "\n")
		estimate <- c(mean(res[-(nrow(res)), i]),
			var(res[-(nrow(res)), i]))
		names(estimate) <- c("mean of simulation",
			"variance of simulation")
		lres[[i]] <- list(statistic=statistic, parameter=parameter,
			method=method, data.name=data.name, p.value=pval, 
			alternative=alternative, estimate=estimate, res=res[,i])
		class(lres[[i]]) <- c("htest", "mc.sim")
		
	}
	class(lres) <- "jclist"
	lres
}


