# Copyright 1998-2002 by Roger Bivand
#

print.sarlm <- function(x, ...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("Type:", x$type, "\n")
	cat("\nCoefficients:\n")
	print(x$coefficients)
	if (x$type == "error") {
		cat("\nLambda:", x$lambda, "\n")
	} else {
		cat("\nRho:", x$rho, "\n")
	}
	cat("\nLog likelihood:", x$LL, "\n")
	invisible(x)
}

summary.sarlm <- function(object, ...)
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
	rownames(object$Coef) <- names(object$coefficients)

	structure(object, class=c("summary.sarlm", class(object)))
}

print.summary.sarlm <- function(x, digits = max(5, .Options$digits - 3),
	signif.stars = FALSE, ...)
{
	cat("\nCall:", paste(deparse(x$call), sep = "", collapse = ""), 
		sep = "", fill=TRUE)
	cat("\nResiduals:\n")
	resid <- x$residuals
	nam <- c("Min", "1Q", "Median", "3Q", "Max")
	rq <- if (length(dim(resid)) == 2) 
		structure(apply(t(resid), 1, quantile), dimnames = list(nam, 
			dimnames(resid)[[2]]))
	else structure(quantile(resid), names = nam)
	print(rq, digits = digits, ...)
	cat("\nType:", x$type, "\n")
	if (x$zero.policy) {
		zero.regs <- attr(ret, "zero.regs")
		if (!is.null(zero.regs))
			cat("Regions with no neighbours included:\n",
			zero.regs, "\n")
	}
	cat("Coefficients:", x$coeftitle, "\n")
	print.coefmat(x$Coef, signif.stars=signif.stars, digits=digits,
		na.print="")
	res <- LR.sarlm(x, x$lm.model)
	if (x$type == "error") {
		cat("\nLambda:", format(signif(x$lambda, digits)),
			"LR test value:", format(signif(res$statistic, digits)),
			"p-value:", format.pval(res$p.value, digits), "\n")
		if (x$ase) cat("Asymptotic standard error:", 
			format(signif(x$lambda.se, digits)),
			"z-value:",format(signif((x$lambda/
				x$lambda.se), digits)),
			"p-value:", format.pval(2*(1-pnorm(abs(x$lambda/
				x$lambda.se))), digits), "\n")
	} else {
		cat("\nRho:", format(signif(x$rho, digits)),
			"LR test value:", format(signif(res$statistic, digits)),
			"p-value:", format.pval(res$p.value, digits), "\n")
		if (x$ase) cat("Asymptotic standard error:", 
			format(signif(x$rho.se, digits)),
			"z-value:",format(signif((x$rho/x$rho.se), digits)),
			"p-value:", format.pval(2 * (1 - pnorm(abs(x$rho/
				x$rho.se))), digits), "\n")
	}
	cat("\nLog likelihood:", x$LL, "for", x$type, "model\n")
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
	cat("\n")
	invisible(x)
}
