# Copyright 2005-8 by Luc Anselin and Roger Bivand
# Kelejian-Prucha generalized moments equations
# for spatial SAR error model
# main function
# Usage:
#    GMerrorsar(formula, data = list(), listw, na.action=na.fail, zero.policy=FALSE, control=list())
# Arguments:
#    formula: standard model formula
#    data: which data frame to search for model variables
#    listw: spatial weights file as list object
#    na.action: standard value
#    zero.policy: allow no-neighbour observations if TRUE
#    control: list of control arguments to optim (such as list(trace=1))
# Details:
#    initializes with ols, calls helper function kpwuwu to build
#    the G and g matrices, calls optim unconstrained optimizer with
#    kpgm as function and plausible starting values to get estimate
#    for lambda, then finds results with spatially weighted least squares
#    and finds LL using Matrix functions
# Value:
# an S3 "gmsar" object

GMerrorsar <- function(#W, y, X, 
	formula, data = list(), listw, na.action=na.fail, 
	zero.policy=NULL, return_LL=FALSE, method="nlminb", 
        control=list(), pars, verbose=NULL, sparse_method="Matrix",
        returnHcov=FALSE, pWOrder=250, tol.Hcov=1.0e-10) {
#	ols <- lm(I(y) ~ I(X) - 1)
        if (is.null(verbose)) verbose <- get("verbose", env = .spdepOptions)
        stopifnot(is.logical(verbose))
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", env = .spdepOptions)
        stopifnot(is.logical(zero.policy))
	mt <- terms(formula, data = data)
	mf <- lm(formula, data, na.action=na.action, method="model.frame")
	na.act <- attr(mf, "na.action")
	if (!is.null(na.act)) {
	    subset <- !(1:length(listw$neighbours) %in% na.act)
	    listw <- subset(listw, subset, zero.policy=zero.policy)
	}

	if (!inherits(listw, "listw")) stop("No neighbourhood list")
	can.sim <- as.logical(NA)
	if (listw$style %in% c("W", "S")) 
		can.sim <- can.be.simmed(listw)

	y <- model.extract(mf, "response")
	if (any(is.na(y))) stop("NAs in dependent variable")
	x <- model.matrix(mt, mf)
	if (any(is.na(x))) stop("NAs in independent variable")
	if (NROW(x) != length(listw$neighbours))
	    stop("Input data and neighbourhood list have different dimensions")

	# added aliased after trying boston with TOWN dummy
	lm.base <- lm(y ~ x - 1)
	aliased <- is.na(coefficients(lm.base))
	cn <- names(aliased)
	names(aliased) <- substr(cn, 2, nchar(cn))
	if (any(aliased)) {
		nacoef <- which(aliased)
		x <- x[,-nacoef]
	}
	ols <- lm(y ~ x - 1)
	if (missing(pars)) {
 	    ubase <- residuals(ols)
	    scorr <- c(crossprod(lag.listw(listw, ubase,
                zero.policy=zero.policy), ubase) / crossprod(ubase, ubase))
            scorr <- scorr / (sum(unlist(listw$weights)) / length(ubase))
            pars <- c(scorr, deviance(ols)/df.residual(ols))
        }
        if (length(pars) !=2 || !is.numeric(pars))
            stop("invalid starting parameter values")
	vv <- .kpwuwu(listw, residuals(ols), zero.policy=zero.policy)
#	nlsres <- nlm(.kpgm, pars, print.level=print.level, gradtol=gradtol, steptol=steptol, iterlim=iterlim, v=vv, verbose=verbose)
#	lambda <- nlsres$estimate[1]
        if (method == "nlminb")
            optres <- nlminb(pars, .kpgm, v=vv, verbose=verbose,
               control=control)
        else 
	    optres <- optim(pars, .kpgm, v=vv, verbose=verbose,
                method=method, control=control)
        if (optres$convergence != 0)
            warning(paste("convergence failure:", optres$message))
	lambda <- optres$par[1]
	names(lambda) <- "lambda"

	wy <- lag.listw(listw, y, zero.policy=zero.policy)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
	n <- NROW(x)
	m <- NCOL(x)
	xcolnames <- colnames(x)
	K <- ifelse(xcolnames[1] == "(Intercept)", 2, 1)
	if (any(is.na(wy)))
	    stop("NAs in lagged dependent variable")
	if (m > 1) {
	    WX <- matrix(nrow=n,ncol=(m-(K-1)))
	    for (k in K:m) {
		wx <- lag.listw(listw, x[,k], zero.policy=zero.policy)
		if (any(is.na(wx)))
		    stop("NAs in lagged independent variable")
		WX[,(k-(K-1))] <- wx
	    }
	}
	if (K == 2) {
# modified to meet other styles, email from Rein Halbersma
		wx1 <- as.double(rep(1, n))
		wx <- lag.listw(listw, wx1, zero.policy=zero.policy)
		if (m > 1) WX <- cbind(wx, WX)
		else WX <- matrix(wx, nrow=n, ncol=1)
	}
	colnames(WX) <- xcolnames
	rm(wx)
	lm.target <- lm(I(y - lambda*wy) ~ I(x - lambda*WX) - 1)
	r <- as.vector(residuals(lm.target))
	fit <- as.vector(y - r)
	p <- lm.target$rank
	SSE <- deviance(lm.target)
	s2 <- SSE/n
	rest.se <- (summary(lm.target)$coefficients[,2])*sqrt((n-p)/n)
	coef.lambda <- coefficients(lm.target)
	names(coef.lambda) <- xcolnames
	call <- match.call()
	names(r) <- names(y)
	names(fit) <- names(y)
	LL <- NULL
	if (return_LL) {
    		if (listw$style %in% c("W", "S") && !can.sim) {
			warning("No log likelihood value available")
		} else {
			if (sparse_method == "spam") {
			  if (listw$style %in% c("W", "S") & can.sim) {
			    csrw <- listw2U_spam(similar.listw_spam(listw))
			  } else csrw <- as.spam.listw(listw)
			  I <- diag.spam(1, n, n)
			} else if (sparse_method == "Matrix") {
			  if (listw$style %in% c("W", "S") & can.sim) {
			    csrw <- listw2U_Matrix(similar.listw_Matrix(listw))
			    similar <- TRUE
			  } else csrw <- as_dsTMatrix_listw(listw)
			  csrw <- as(csrw, "CsparseMatrix")
			  I <- as_dsCMatrix_I(n)
			} else stop("unknown sparse_method")
			gc(FALSE)
			yl <- y - lambda*wy
			xl <- x - lambda*WX
			xl.q <- qr.Q(qr(xl))
			xl.q.yl <- t(xl.q) %*% yl
			SSE <- t(yl) %*% yl - t(xl.q.yl) %*% xl.q.yl
			s2 <- SSE/n
			if (sparse_method == "spam") {
			  Jacobian <- determinant((I - lambda * csrw), 
			    logarithm=TRUE)$modulus
			} else if (sparse_method == "Matrix") {
                             .f <- if (package_version(packageDescription(
                                 "Matrix")$Version) > "0.999375-30") 2 else 1
			  Jacobian <- .f * determinant(I - lambda * csrw,
 			    logarithm=TRUE)$modulus
			}
			gc(FALSE)
			LL <- (Jacobian -
				((n/2)*log(2*pi)) - (n/2)*log(s2) - 
				(1/(2*(s2)))*SSE)
		}
	}
        Hcov <- NULL
        if (returnHcov) {
            W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
            pp <- ols$rank
            p1 <- 1L:pp
            R <- chol2inv(ols$qr$qr[p1, p1, drop = FALSE])
            B <- tcrossprod(R, x)
            B <- as(powerWeights(W=W, rho=lambda, order=pWOrder,
                X=B, tol=tol.Hcov), "matrix")
            C <- x %*% R
            C <- as(powerWeights(W=t(W), rho=lambda, order=pWOrder,
                X=C, tol=tol.Hcov), "matrix")
            Hcov <- B %*% C
            attr(Hcov, "method") <- "Matrix"
        }

	ret <- structure(list(lambda=lambda,
		coefficients=coef.lambda, rest.se=rest.se, 
		s2=s2, SSE=SSE, parameters=(m+2), lm.model=ols, 
		call=call, residuals=r, lm.target=lm.target,
		fitted.values=fit, formula=formula, aliased=aliased,
		zero.policy=zero.policy, LL=LL, vv=vv, optres=optres,
                pars=pars, Hcov=Hcov), class=c("gmsar"))
	if (zero.policy) {
		zero.regs <- attr(listw$neighbours, 
			"region.id")[which(card(listw$neighbours) == 0)]
		if (length(zero.regs) > 0)
			attr(ret, "zero.regs") <- zero.regs
	}
	if (!is.null(na.act))
		ret$na.action <- na.act
	ret
}

# Copyright 2005 by Roger Bivand

residuals.gmsar <- function(object, ...) {
	if (is.null(object$na.action))
		object$residuals
	else napredict(object$na.action, object$residuals)
}

deviance.gmsar <- function(object, ...) {
	deviance(object$lm.target)
}


coef.gmsar <- function(object, ...) {
	ret <- c(object$coefficients, object$lambda)
	ret
}

fitted.gmsar <- function(object, ...) {
	if (is.null(object$na.action))
		object$fitted.values
	else napredict(object$na.action, object$fitted.values)
}


print.gmsar <- function(x, ...)
{
	cat("\nCall:\n")
	print(x$call)
	cat("\n")
	cat("\nCoefficients:\n")
	print(coef(x))
	cat("\nLog likelihood:", logLik(x), "\n")
	invisible(x)
}

summary.gmsar <- function(object, correlation = FALSE, Hausman=FALSE, ...)
{
	object$coeftitle <- "(GM standard errors)"
	object$Coef <- cbind(object$coefficients, object$rest.se, 
		object$coefficients/object$rest.se,
		2*(1-pnorm(abs(object$coefficients/object$rest.se))))
	colnames(object$Coef) <- c("Estimate", "Std. Error", 
		"z value", "Pr(>|z|)")
	if (is.null(object$LL)) object$LR1 <- NULL
	else object$LR1 <- LR1.gmsar(object)
	rownames(object$Coef) <- names(object$coefficients)
        if (Hausman && !is.null(object$Hcov)) {
                object$Haus <- Hausman.test(object)
        }

	structure(object, class=c("summary.gmsar", class(object)))
}

LR1.gmsar <- function(object)
{
	if (!inherits(object, "gmsar")) stop("Not a gmsar object")
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
	attr(estimate, "names") <- c("Log likelihood of GM estimator", "Log likelihood of OLS fit")
	method <- "Likelihood Ratio diagnostics for spatial dependence"
	res <- list(statistic=statistic, parameter=parameter,
		p.value=p.value, estimate=estimate, method=method)
	class(res) <- "htest"
	res
}

logLik.gmsar <- function(object, ...) {
	if (is.null(object$LL)) {
		warning("Model fitted without LL")
		return(NULL)
	}
	LL <- c(object$LL)
	class(LL) <- "logLik"
	N <- length(residuals(object))
	attr(LL, "nall") <- N
	attr(LL, "nobs") <- N
	attr(LL, "df") <- object$parameters
	LL
}

print.summary.gmsar <- function(x, digits = max(5, .Options$digits - 3),
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
	cat("\nType: GM SAR estimator\n")
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

	res <- x$LR1
	cat("\nLambda:", format(signif(x$lambda, digits)))
	if (!is.null(res)) cat(" LR test value:", format(signif(res$statistic, 
		digits)), "p-value:", format.pval(res$p.value, digits), "\n")
	else cat("\n")
	if (!is.null(x$LL)) 
		cat("\nLog likelihood:", logLik(x), "for GM model\n")
	cat("ML residual variance (sigma squared): ", 
		format(signif(x$s2, digits)), ", (sigma: ", 
		format(signif(sqrt(x$s2), digits)), ")\n", sep="")
	cat("Number of observations:", length(x$residuals), "\n")
	cat("Number of parameters estimated:", x$parameters, "\n")
	if (!is.null(res)) cat("AIC: ", format(signif(AIC(x), digits)), 
		", (AIC for lm: ", format(signif(AIC(x$lm.model), digits)),
		")\n", sep="")
	if (!is.null(x$Haus)) {
	    cat("Hausman test: ", format(signif(x$Haus$statistic, 
		digits)), ", df: ", format(x$Haus$parameter),
                       ", p-value: ", format.pval(x$Haus$p.value, digits),
                       "\n", sep="")
	}
    	cat("\n")
        invisible(x)
}

# Copyright 2004 by Luc Anselin
# Kelejian-Prucha generalized moments equations
# helper function to provide function to nonlinear optimizer
# must have parameter vector first for nlm
# Usage:
#    kpgm(par,v)
# Arguments:
#    par: 2x1 parameter vector rho,sig2
#    v: list containing bigG and litg as computed by kpwuwu
# Details:
#    sets up the equation as squared residuals
# Value:
#    value: evaluated nonlinear least squares for parameter value

.kpgm <- function(rhopar,v,verbose=FALSE) {
  vv <- v$bigG %*% c(rhopar[1],rhopar[1]^2,rhopar[2]) - v$litg
  value <- sum(vv^2)
  if (verbose)
    cat("function:", value, "lambda:", rhopar[1], "sig2:", rhopar[2], "\n")
  value
  
}


# Copyright 2004 by Luc Anselin
# Kelejian-Prucha generalized moments equations
# helper function
# Usage:
#    kpwuwu(listw,u)
# Arguments:
#    listw: spatial weights file as listw object
#    u: OLS residual vector
#    zero.policy: allow no-neighbour observations if TRUE
# Details:
#    sets up the bigG matrix and littleg vector needed
#    for the nonlinear least squares in the GM estimator
#    see Kelejian-Prucha(1999) p. 515
# Value:
# a list with two elements
#    bigG: the 3x3 G matrix
#    litg: the 3x1 g vector

.kpwuwu <- function(W, u, zero.policy=FALSE) {
	n <- length(u)
# Gianfranco Piras 081119 
        trwpw <- sum(unlist(W$weights)^2)
#	tt <- matrix(0,n,1)
#	for (i in 1:n) {tt[i] <- sum(W$weights[[i]]^2) }
#	trwpw <- sum(tt)
	wu <- lag.listw(W, u, zero.policy=zero.policy)
	wwu <- lag.listw(W, wu, zero.policy=zero.policy)
    	uu <- crossprod(u,u)
    	uwu <- crossprod(u,wu)
    	uwpuw <- crossprod(wu,wu)
    	uwwu <- crossprod(u,wwu)
    	wwupwu <- crossprod(wwu,wu)
    	wwupwwu <- crossprod(wwu,wwu)
    	bigG <- matrix(0,3,3)
    	bigG[,1] <- c(2*uwu,2*wwupwu,(uwwu+uwpuw))/n
    	bigG[,2] <- - c(uwpuw,wwupwwu,wwupwu) / n
    	bigG[,3] <- c(1,trwpw/n,0)
    	litg <- c(uu,uwpuw,uwu) / n
    	list(bigG=bigG,litg=litg)
}

