# Copyright 1998-2004 by Roger Bivand (Wald test suggested by Rein Halbersma,
# output of correlations suggested by Michael Tiefelsdorf)
#

print.sarlm <- function(x, ...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("Type:", x$type, "\n")
	cat("\nCoefficients:\n")
	print(coef(x))
	cat("\nLog likelihood:", logLik(x), "\n")
	invisible(x)
}

summary.sarlm <- function(object, correlation = FALSE, ...)
{
	if (object$type == "error" || ((object$type == "lag" || 
		object$type == "mixed") && object$ase)) {
		object$coeftitle <- "(asymptotic standard errors)"
		object$Coef <- cbind(object$coefficients, object$rest.se, 
			object$coefficients/object$rest.se,
			2*(1-pnorm(abs(object$coefficients/object$rest.se))))
		colnames(object$Coef) <- c("Estimate", "Std. Error", 
			"z value", "Pr(>|z|)")
	} else {
	    # intercept-only bug fix Larry Layne 20060404
	    if (is.null(object$LLs)) {
		object$Coef <- cbind(object$coefficients)
		colnames(object$Coef) <- c("Estimate")

	    } else {
		object$coeftitle <- "(log likelihood/likelihood ratio)"
		m <- length(object$coefficients)
		LLs <- numeric(m)
		LRs <- numeric(m)
		Pvals <- numeric(m)
		if (length(object$LLs) == m) {
			for (i in 1:m) {
				LLs[i] <- c(object$LLs[[i]])
				res <- LR.sarlm(object, object$LLs[[i]])
				LRs[i] <- res$statistic
				Pvals[i] <- res$p.value
			}
		} else {
			LLs[1] <- LRs[1] <- Pvals[1] <- NA
			for (i in 1:(m-1)) {
				LLs[i+1] <- c(object$LLs[[i]])
				res <- LR.sarlm(object, object$LLs[[i]])
				LRs[i+1] <- res$statistic
				Pvals[i+1] <- res$p.value
			}
		}
		object$Coef <- cbind(object$coefficients, LLs, LRs, Pvals)
		colnames(object$Coef) <- c("Estimate", "Log likelihood",
			"LR statistic", "Pr(>|z|)")
	    }
	}
	if (object$ase) {
		object$Wald1 <- Wald1.sarlm(object)
		if (correlation) {
			object$correlation <- diag((diag(object$resvar))
				^(-1/2)) %*% object$resvar %*% 
				diag((diag(object$resvar))^(-1/2))
			dimnames(object$correlation) <- dimnames(object$resvar)
		}
	}
	object$LR1 <- LR1.sarlm(object)
	rownames(object$Coef) <- names(object$coefficients)

	structure(object, class=c("summary.sarlm", class(object)))
}

LR1.sarlm <- function(object)
{
	if (!inherits(object, "sarlm")) stop("Not a sarlm object")
	LLx <- logLik(object)
	LLy <- logLik(object$lm.model)
	statistic <- 2*(LLx - LLy)
	attr(statistic, "names") <- "Likelihood ratio"
	parameter <- abs(attr(LLx, "df") - attr(LLy, "df"))
	if (parameter < 1) 
		stop("non-positive degrees of freedom: no test possible")
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	estimate <- c(LLx, LLy)
	if (object$type == "error") alt <- "spatial error model"
	else alt <- "spatial lag model"
	attr(estimate, "names") <- c(paste("Log likelihood of",
		alt), paste("Log likelihood of OLS fit",
		deparse(substitute(y))))
	method <- "Likelihood Ratio diagnostics for spatial dependence"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res
}

Wald1.sarlm <- function(object) {
	if (!inherits(object, "sarlm")) stop("Not a sarlm object")
	if (!object$ase) 
		stop("Cannot compute Wald statistic: parameter a.s.e. missing")
	LLx <- logLik(object)
	LLy <- logLik(object$lm.model)
	if (object$type == "lag" || object$type == "mixed") {
		estimate <- object$rho
		statistic <- (object$rho / object$rho.se)^2
	} else {
		estimate <- object$lambda
		statistic <- (object$lambda / object$lambda.se)^2
	}
	attr(statistic, "names") <- "Wald statistic"
	parameter <- abs(attr(LLx, "df") - attr(LLy, "df"))
	if (parameter < 1) 
		stop("non-positive degrees of freedom: no test possible")
	attr(parameter, "names") <- "df"
	p.value <- 1 - pchisq(abs(statistic), parameter)
	method <- "Wald diagnostics for spatial dependence"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res

}


print.summary.sarlm <- function(x, digits = max(5, .Options$digits - 3),
	signif.stars = FALSE, ...)
{
	cat("\nCall:", deparse(x$call),	sep = "", fill=TRUE)
	cat("\nResiduals:\n")
	resid <- residuals(x)
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2) 
		structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
			dimnames(resid)[[2]]))
	else structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
	cat("\nType:", x$type, "\n")
	if (x$zero.policy) {
		zero.regs <- attr(x, "zero.regs")
		if (!is.null(zero.regs))
			cat("Regions with no neighbours included:\n",
			zero.regs, "\n")
	}
	cat("Coefficients:", x$coeftitle, "\n")
	coefs <- x$Coef
	if (!is.null(aliased <- x$aliased) && any(x$aliased)){
		cat("    (", table(aliased)["TRUE"], 
			" not defined because of singularities)\n", sep = "")
		cn <- names(aliased)
		coefs <- matrix(NA, length(aliased), 4, dimnames = list(cn, 
                	colnames(x$Coef)))
            	coefs[!aliased, ] <- x$Coef
	}
	printCoefmat(coefs, signif.stars=signif.stars, digits=digits,
		na.print="NA")
#	res <- LR.sarlm(x, x$lm.model)
	res <- x$LR1
	if (x$type == "error") {
		cat("\nLambda:", format(signif(x$lambda, digits)),
			"LR test value:", format(signif(res$statistic, digits)),
			"p-value:", format.pval(res$p.value, digits), "\n")
		if (x$ase) {
			cat("Asymptotic standard error:", 
			format(signif(x$lambda.se, digits)),
			"z-value:",format(signif((x$lambda/
				x$lambda.se), digits)),
			"p-value:", format.pval(2*(1-pnorm(abs(x$lambda/
				x$lambda.se))), digits), "\n")
			cat("Wald statistic:", format(signif(x$Wald1$statistic, 
			digits)), "p-value:", format.pval(x$Wald1$p.value, 
			digits), "\n")
		}
	} else {
		cat("\nRho:", format(signif(x$rho, digits)),
			"LR test value:", format(signif(res$statistic, digits)),
			"p-value:", format.pval(res$p.value, digits), "\n")
		if (x$ase) {
			cat("Asymptotic standard error:", 
			format(signif(x$rho.se, digits)),
			"z-value:",format(signif((x$rho/x$rho.se), digits)),
			"p-value:", format.pval(2 * (1 - pnorm(abs(x$rho/
				x$rho.se))), digits), "\n")
			cat("Wald statistic:", format(signif(x$Wald1$statistic, 
			digits)), "p-value:", format.pval(x$Wald1$p.value, 
			digits), "\n")
		}

	}
	cat("\nLog likelihood:", logLik(x), "for", x$type, "model\n")
	cat("ML residual variance (sigma squared): ", 
		format(signif(x$s2, digits)), ", (sigma: ", 
		format(signif(sqrt(x$s2), digits)), ")\n", sep="")
	cat("Number of observations:", length(x$residuals), "\n")
	cat("Number of parameters estimated:", x$parameters, "\n")
	cat("AIC: ", format(signif(AIC(x), digits)), ", (AIC for lm: ",
		format(signif(AIC(x$lm.model), digits)), ")\n", sep="")
	if (x$type != "error" && x$ase) {
		cat("LM test for residual autocorrelation\n")
		cat("test value:", format(signif(x$LMtest, digits)),
			"p-value:", format.pval((1 - pchisq(x$LMtest, 1)), 
			digits), "\n")
	}
    	correl <- x$correlation
    	if (!is.null(correl)) {
        	p <- NCOL(correl)
        	if (p > 1) {
            		cat("\nCorrelation of Coefficients:\n")
                	correl <- format(round(correl, 2), nsmall = 2, 
                  	digits = digits)
                	correl[!lower.tri(correl)] <- ""
                	print(correl[-1, -p, drop = FALSE], quote = FALSE)
            	}
    	}
    	cat("\n")
        invisible(x)
}
