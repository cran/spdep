# Copyright 2001-5 by Roger Bivand 
#

moran <- function(x, listw, n, S0, zero.policy=FALSE, NAOK=FALSE) {
	n1 <- length(listw$neighbours)
	x <- c(x)
	if (n1 != length(x)) stop("objects of different length")
	xx <- mean(x, na.rm=NAOK)
	z <- x - xx
	zz <- sum(z^2, na.rm=NAOK)
	K <- (length(x)*sum(z^4, na.rm=NAOK))/(zz^2)
	lz <- lag.listw(listw, z, zero.policy=zero.policy, NAOK=NAOK)
#	I <- (n / S0) * ((t(z) %*% lz) / zz)
	I <- (n / S0) * ((sum(z*lz, na.rm=NAOK)) / zz)
	res <- list(I=I, K=K)
	res
}

moran.test <- function(x, listw, randomisation=TRUE, zero.policy=FALSE,
	alternative="greater", rank = FALSE, na.action=na.fail, spChk=NULL) {
	alternative <- match.arg(alternative, c("greater", "less", "two.sided"))
	if (!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
#	if (any(is.na(x))) stop("NA in X")
	xname <- deparse(substitute(x))
	wname <- deparse(substitute(listw))
	NAOK <- deparse(substitute(na.action)) == "na.pass"
	x <- na.action(x)
	na.act <- attr(x, "na.action")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	
	wc <- spweights.constants(listw, zero.policy=zero.policy)
	S02 <- wc$S0*wc$S0
	res <- moran(x, listw, wc$n, wc$S0, zero.policy=zero.policy, 
		NAOK=NAOK)
	I <- res$I
	K <- res$K
	if (rank) K <- (3*(3*wc$n^2 -7))/(5*(wc$n^2 - 1))
	EI <- (-1) / wc$n1
	if(randomisation) {
		VI <- wc$n*(wc$S1*(wc$nn - 3*wc$n + 3) - wc$n*wc$S2 + 3*S02)
		tmp <- K*(wc$S1*(wc$nn - wc$n) - 2*wc$n*wc$S2 + 6*S02)
		VI <- (VI - tmp) / (wc$n1*wc$n2*wc$n3*S02)
		VI <- VI - EI^2
	} else {
		VI <- (wc$nn*wc$S1 - wc$n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
		VI <- VI - EI^2
	}
	ZI <- (I - EI) / sqrt(VI)
	statistic <- ZI
	names(statistic) <- "Moran I statistic standard deviate"
        if (alternative == "two.sided") 
		PrI <- 2 * pnorm(abs(ZI), lower.tail=FALSE)
        else if (alternative == "greater")
            PrI <- pnorm(ZI, lower.tail=FALSE)
        else PrI <- pnorm(ZI)
	if (PrI < 0 || PrI > 1) 
		warning("Out-of-range p-value: reconsider test arguments")
	vec <- c(I, EI, VI)
	names(vec) <- c("Moran I statistic", "Expectation", "Variance")
	method <- paste("Moran's I test under", ifelse(randomisation,
	    "randomisation", "normality"))
	data.name <- paste(xname, ifelse(rank,
		"using rank correction",""), "\nweights:",
		wname, ifelse(is.null(na.act), "", paste("\nomitted:", 
	    paste(na.act, collapse=", "))), "\n")
	res <- list(statistic=statistic, p.value=PrI, estimate=vec, 
	    alternative=alternative, method=method, data.name=data.name)
	if (!is.null(na.act)) attr(res, "na.action") <- na.act
	class(res) <- "htest"
	res
}

moran.mc <- function(x, listw, nsim, zero.policy=FALSE,
	alternative="greater", na.action=na.fail, spChk=NULL) {
	alternative <- match.arg(alternative, c("greater", "less"))
	if(!inherits(listw, "listw")) stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if(missing(nsim)) stop("nsim must be given")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	cards <- card(listw$neighbours)
	if (!zero.policy && any(cards == 0))
		stop("regions with no neighbours found")
#	if (any(is.na(x))) stop("NA in X")
	xname <- deparse(substitute(x))
	wname <- deparse(substitute(listw))
	if (deparse(substitute(na.action)) == "na.pass")
	    stop("na.pass not permitted")
	x <- na.action(x)
	na.act <- attr(x, "na.action")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	if (nsim > gamma(n+1)) stop("nsim too large for this n")
	if (nsim < 1) stop("nsim too small")
	
	S0 <- Szero(listw)
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) res[i] <- moran(sample(x), listw, n, S0,
	    zero.policy)$I
	res[nsim+1] <- moran(x, listw, n, S0, zero.policy)$I
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
	method <- "Monte-Carlo simulation of Moran's I"
	data.name <- paste(xname, "\nweights:",
	    wname, ifelse(is.null(na.act), "", paste("\nomitted:", 
	    paste(na.act, collapse=", "))), "\nnumber of simulations + 1:",
	    nsim+1, "\n")
	lres <- list(statistic=statistic, parameter=parameter,
	    p.value=pval, alternative=alternative, method=method, 
	    data.name=data.name, res=res)
	if (!is.null(na.act)) attr(lres, "na.action") <- na.act
	class(lres) <- c("htest", "mc.sim")
	lres
}


