# Copyright 1998-2009 by Roger Bivand (Wald test suggested by Rein Halbersma,
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

summary.sarlm <- function(object, correlation = FALSE, Nagelkerke=FALSE,
 Hausman=FALSE, ...)
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
            if (!is.null(object$rest.se)) {
		object$coeftitle <- "(numerical Hessian approximate standard errors)"
		object$Coef <- cbind(object$coefficients, object$rest.se, 
			object$coefficients/object$rest.se,
			2*(1-pnorm(abs(object$coefficients/object$rest.se))))
		colnames(object$Coef) <- c("Estimate", "Std. Error", 
			"z value", "Pr(>|z|)")
	        rownames(object$Coef) <- names(object$coefficients)
              }
	}

        if (Nagelkerke) {
            nk <- NK.sarlm(object)
            if (!is.null(nk)) object$NK <- nk
        }
        if (Hausman && object$type == "error" && !is.null(object$Hcov)) {
                object$Haus <- Hausman.test(object)
        }
	if (object$type == "error") {
		object$Wald1 <- Wald1.sarlm(object)
		if (correlation) {
                        oresvar <- object$resvar
                        ctext <- "Correlation of coefficients"
                        if (is.null(oresvar) || is.logical(oresvar)) {
                            oresvar <- object$fdHess
                            ctext <- ifelse(object$insert,
                                "Approximate correlation of coefficients",
                                "** Guesswork correlation of coefficients **")
                        }
			object$correlation <- diag((diag(oresvar))
				^(-1/2)) %*% oresvar %*% 
				diag((diag(oresvar))^(-1/2))
			dimnames(object$correlation) <- dimnames(oresvar)
                        object$correltext <- ctext
		}
	} else if (object$type != "error") {
		object$Wald1 <- Wald1.sarlm(object)
		if (correlation) {
                        oresvar <- object$resvar
                        ctext <- "Correlation of coefficients"
                        if (is.null(oresvar) || is.logical(oresvar)) {
                            oresvar <- object$fdHess
                            ctext <- "Approximate correlation of coefficients"
                        }
			object$correlation <- diag((diag(oresvar))
				^(-1/2)) %*% oresvar %*% 
				diag((diag(oresvar))^(-1/2))
			dimnames(object$correlation) <- dimnames(oresvar)
                        object$correltext <- ctext
		}
        }
	object$LR1 <- LR1.sarlm(object)

	structure(object, class=c("summary.sarlm", class(object)))
}

NK.sarlm <- function(obj) {
     n <- length(obj$residuals)
     nullLL <- obj$LLNullLlm
     if (is.null(nullLL)) return(nullLL)
     c(1 - exp(-(2/n)*(logLik(obj) - nullLL)))
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
#	if (!object$ase) 
#		stop("Cannot compute Wald statistic: parameter a.s.e. missing")
	LLx <- logLik(object)
	LLy <- logLik(object$lm.model)
	if (object$type == "lag" || object$type == "mixed") {
		estimate <- object$rho
                rse <- object$rho.se
                if (is.null(rse)) return(rse)
		statistic <- (object$rho / rse)^2
		attr(statistic, "names") <- ifelse(is.logical(fdHess), 
                    "Wald statistic", "Approximate Wald statistic")
	} else {
		estimate <- object$lambda
                lse <- object$lambda.se
                if (is.null(lse)) return(lse)
		statistic <- (object$lambda / lse)^2
		attr(statistic, "names") <- ifelse(is.logical(fdHess), 
                    "Wald statistic", "Approximate Wald statistic")
	}
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

Hausman.test <- function(object, ...)
    UseMethod("Hausman.test", object)

Hausman.test.sarlm <- function(object, ..., tol=NULL) {
    if (!inherits(object, "sarlm")) stop("not a sarlm object")
    if (object$type != "error") stop("not a spatial error model")
    fmeth <- ifelse(object$method != "eigen", "(approximate)", "(asymptotic)") 
    if (is.null(object$Hcov)) stop("Vo not available")
    s2 <- object$s2
    Vo <- s2 * object$Hcov
    Vs <- s2 * summary.lm(object$lm.target, corr = FALSE)$cov.unscaled
    d <- coef(object$lm.model) - coef(object$lm.target)
    if (!is.null(tol)) VV <- try(solve((Vo - Vs), tol=tol))
    else VV <- try(solve(Vo - Vs))
    if (class(VV) == "try.error") {
        warning("(Vo - Vs) inversion failure")
        return(NULL)
    }
    statistic <- t(d) %*% VV %*% d
    attr(statistic, "names") <- "Hausman test"
    parameter <- length(d)
    attr(parameter, "names") <- "df"
    p.value <- 1 - pchisq(abs(statistic), parameter)
    method <- paste("Spatial Hausman test", fmeth)
    data.name <- strwrap(deparse(object$formula), exdent=4)
    if (length(data.name) > 1) 
        data.name <- paste(data.name, collapse="\n    ")
    res <- list(statistic = statistic, parameter = parameter, 
        p.value = p.value, method = method, data.name=data.name)
    class(res) <- "htest"
    res
}

Hausman.test.gmsar <- function(object, ..., tol=NULL) {
    if (!inherits(object, "gmsar")) stop("not a gmsar object")
    if (is.null(object$Hcov)) stop("Vo not available")
    fmeth <- "(approximate)"
    s2 <- object$s2
    Vo <- s2 * object$Hcov
    Vs <- s2 * summary.lm(object$lm.target, corr = FALSE)$cov.unscaled
    d <- coef(object$lm.model) - coef(object$lm.target)
    if (!is.null(tol)) VV <- try(solve((Vo - Vs), tol=tol))
    else VV <- try(solve(Vo - Vs))
    if (class(VV) == "try.error") {
        warning("(Vo - Vs) inversion failure")
        return(NULL)
    }
    statistic <- t(d) %*% VV %*% d
    attr(statistic, "names") <- "Hausman test"
    parameter <- length(d)
    attr(parameter, "names") <- "df"
    p.value <- 1 - pchisq(abs(statistic), parameter)
    method <- paste("Spatial Hausman test", fmeth)
    data.name <- strwrap(deparse(object$formula), exdent=4)
    if (length(data.name) > 1) 
        data.name <- paste(data.name, collapse="\n    ")
    res <- list(statistic = statistic, parameter = parameter, 
        p.value = p.value, method = method, data.name=data.name)
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
        if (!is.null(x$coeftitle)) {
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
	}
#	res <- LR.sarlm(x, x$lm.model)
	res <- x$LR1
	if (x$type == "error") {
		cat("\nLambda: ", format(signif(x$lambda, digits)),
			", LR test value: ", format(signif(res$statistic,
                        digits)), ", p-value: ", format.pval(res$p.value,
                        digits), "\n", sep="")
		if (!is.null(x$lambda.se)) {
                    pref <- ifelse(x$ase, "Asymptotic",
                        "Approximate (numerical Hessian)")
		    cat(pref, " standard error: ", 
		        format(signif(x$lambda.se, digits)),
			"\n    z-value: ",format(signif((x$lambda/
				x$lambda.se), digits)),
			", p-value: ", format.pval(2*(1-pnorm(abs(x$lambda/
				x$lambda.se))), digits), "\n", sep="")
		    cat("Wald statistic: ", format(signif(x$Wald1$statistic, 
			digits)), ", p-value: ", format.pval(x$Wald1$p.value, 
			digits), "\n", sep="")
		}
	} else {
		cat("\nRho: ", format(signif(x$rho, digits)), 
                    ", LR test value: ", format(signif(res$statistic, digits)),
		    ", p-value: ", format.pval(res$p.value, digits), "\n",
                    sep="")
                if (!is.null(x$rho.se)) {
                  pref <- ifelse(x$ase, "Asymptotic",
                    "Approximate (numerical Hessian)")
		  cat(pref, " standard error: ", 
			format(signif(x$rho.se, digits)), "\n    z-value: ", 
			format(signif((x$rho/x$rho.se), digits)),
			", p-value: ", format.pval(2 * (1 - pnorm(abs(x$rho/
				x$rho.se))), digits), "\n", sep="")
                }
		if (!is.null(x$Wald1)) {
		    cat("Wald statistic: ", format(signif(x$Wald1$statistic, 
			digits)), ", p-value: ", format.pval(x$Wald1$p.value, 
			digits), "\n", sep="")
		}

	}
	cat("\nLog likelihood:", logLik(x), "for", x$type, "model\n")
	cat("ML residual variance (sigma squared): ", 
		format(signif(x$s2, digits)), ", (sigma: ", 
		format(signif(sqrt(x$s2), digits)), ")\n", sep="")
        if (!is.null(x$NK)) cat("Nagelkerke pseudo-R-squared:",
            format(signif(x$NK, digits)), "\n")
	cat("Number of observations:", length(x$residuals), "\n")
	cat("Number of parameters estimated:", x$parameters, "\n")
	cat("AIC: ", format(signif(AIC(x), digits)), ", (AIC for lm: ",
		format(signif(AIC(x$lm.model), digits)), ")\n", sep="")
	if (x$type == "error") {
		if (!is.null(x$Haus)) {
		    cat("Hausman test: ", format(signif(x$Haus$statistic, 
			digits)), ", df: ", format(x$Haus$parameter),
                        ", p-value: ", format.pval(x$Haus$p.value, digits),
                        "\n", sep="")
		}
        }
	if (x$type != "error" && x$ase) {
		cat("LM test for residual autocorrelation\n")
		cat("test value: ", format(signif(x$LMtest, digits)),
			", p-value: ", format.pval((1 - pchisq(x$LMtest, 1)), 
			digits), "\n", sep="")
	}
        if (x$type != "error" && !is.null(x$LLCoef)) {
		cat("\nCoefficients: (log likelihood/likelihood ratio)\n")
		printCoefmat(x$LLCoef, signif.stars=signif.stars,
			digits=digits, na.print="NA")
        }
    	correl <- x$correlation
    	if (!is.null(correl)) {
        	p <- NCOL(correl)
        	if (p > 1) {
            		cat("\n", x$correltext, "\n")
                	correl <- format(round(correl, 2), nsmall = 2, 
                  	digits = digits)
                	correl[!lower.tri(correl)] <- ""
                	print(correl[-1, -p, drop = FALSE], quote = FALSE)
            	}
    	}
    	cat("\n")
        invisible(x)
}
