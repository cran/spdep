# Copyright 2002 by Roger Bivand 
#

LR.sarlm <- function(x, y)
{
	if (!inherits(x, "logLik")) LLx <- logLik(x)
	else LLx <- c(x)
	if (!inherits(y, "logLik")) LLy <- logLik(y)
	else LLy <- c(y)
	statistic <- 2*(LLx - LLy)
	attr(statistic, "names") <- "Likelihood ratio"
	parameter <- 1
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	estimate <- c(LLx, LLy)
	attr(estimate, "names") <- c(paste("Log likelihood of",
		deparse(substitute(x))), paste("Log likelihood of",
		deparse(substitute(y))))
	method <- "Likelihood ratio for spatial linear models"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res
}

logLik.sarlm <- function(object, ...) {
	LL <- c(object$LL)
	class(LL) <- "logLik"
	N <- length(residuals(object))
	attr(LL, "nall") <- N
	attr(LL, "nobs") <- N
	attr(LL, "df") <- object$parameters
	LL
}

