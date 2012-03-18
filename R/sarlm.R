# Copyright 2002-12 by Roger Bivand
#

residuals.sarlm <- function(object, ...) {
	if (is.null(object$na.action))
		object$residuals
	else napredict(object$na.action, object$residuals)
}

deviance.sarlm <- function(object, ...) {
	object$SSE
}

coef.sarlm <- function(object, ...) {
	ret <- NULL
#	ret <- sqrt(object$s2)
#	names(ret) <- "sigma"
	if (object$type == "error") ret <- c(ret, object$lambda)
	else if (object$type == "lag" || object$type == "mixed")
            ret <- c(ret, object$rho)
        else if (object$type == "sac" || object$type == "sacmixed")
            ret <- c(ret, object$rho, object$lambda)
	ret <- c(ret, object$coefficients)

	ret
}

vcov.sarlm <- function(object, ...) {
	if (object$ase) res <- object$resvar[-1,-1]
        else {
            if (!is.null(object$fdHess)) {
                if (object$insert) res <- object$resvar[-1,-1]
                else res <- object$resvar
            } else {
                stop("vcov not available for this model")
            }
        }
        res
}


fitted.sarlm <- function(object, ...) {
	if (is.null(object$na.action))
		object$fitted.values
	else napredict(object$na.action, object$fitted.values)
}

predict.sarlm <- function(object, newdata=NULL, listw=NULL, 
	zero.policy=NULL, legacy=TRUE, power=NULL, order=250, tol=.Machine$double.eps^(3/5), pred.se=FALSE, lagImpact=NULL, ...) {
        if (is.null(zero.policy))
            zero.policy <- get("zeroPolicy", envir = .spdepOptions)
        stopifnot(is.logical(zero.policy))
        if (object$type == "sac") stop("no predict method for sac")
        if (is.null(power)) power <- object$method != "eigen"
        stopifnot(is.logical(legacy))
        stopifnot(is.logical(power))
        if (pred.se && object$type == "error") {
            pred.se <- FALSE
            warning("standard error estimates not available for error models")
        }
        if (pred.se && is.null(lagImpact))
            stop("lagImpact object from impact method required for standard error estimate")
	if (is.null(newdata)) {
		res <- fitted.values(object)
		X <- object$X
		B <- object$coefficients
		y <- object$y
		tarX <- object$tarX
		tary <- object$tary
		if (object$type == "error") {
			attr(res, "trend") <- as.vector(X %*% B)
			attr(res, "signal") <- as.vector( -1 * (tary - y) - 					-1 * (tarX - X) %*% B)
		} else {
			attr(res, "trend") <- as.vector(X %*% B)
			attr(res, "signal") <- as.vector( -1 * (tary - y))
		}
	}
	else {
		if (object$type == "error") {
                    if (object$etype == "error") {
			B <- object$coefficients
#			tt <- terms(object$lm.model) 
#			X <- model.matrix(delete.response(tt), data=newdata)
                        frm <- formula(object$call)
			mt <- delete.response(terms(frm, data=newdata))
#			mf <- lm(object$formula, newdata, method="model.frame")
			mf <- model.frame(mt, newdata)
			X <- model.matrix(mt, mf)
#  accommodate aliased coefficients 120314
                        if (any(object$aliased))
                            X <- X[,-which(object$aliased)]
			trend <- X %*% B
			signal <- rep(0, length(trend))
			res <- trend + signal
			attr(res, "trend") <- trend
			attr(res, "signal") <- signal
                    } else if (object$etype == "emixed") {
			if (is.null(listw) || !inherits(listw, "listw")) 
				stop ("spatial weights list required")
			if (nrow(newdata) != length(listw$neighbours))
				stop("mismatch between newdata and spatial weights")
			B <- object$coefficients
#			mt <- terms(object$formula, data = newdata)
                        frm <- formula(object$call)
			mt <- delete.response(terms(frm, data=newdata))
#			mf <- lm(object$formula, newdata, method="model.frame")
			mf <- model.frame(mt, newdata)
			X <- model.matrix(mt, mf)
			K <- ifelse(colnames(X)[1] == "(Intercept)", 2, 1)
			m <- ncol(X)
			WX <- matrix(nrow=nrow(X),ncol=(m-(K-1)))
			for (k in K:m) {
				wx <- lag.listw(listw, X[,k], 
				    zero.policy=zero.policy)
				if (any(is.na(wx))) 
				    stop("NAs in lagged independent variable")
				WX[,(k-(K-1))] <- wx
			}
			X <- cbind(X, WX)
#  accommodate aliased coefficients 120314
                        if (any(object$aliased))
                            X <- X[,-which(object$aliased)]
			trend <- X %*% B
			signal <- rep(0, length(trend))
			res <- trend + signal
			attr(res, "trend") <- trend
			attr(res, "signal") <- signal
                    } else stop("unkown error model etype")
		} else if (object$type == "mixed") {
			if (is.null(listw) || !inherits(listw, "listw")) 
				stop ("spatial weights list required")
			if (nrow(newdata) != length(listw$neighbours))
				stop("mismatch between newdata and spatial weights")
			B <- object$coefficients
#			mt <- terms(object$formula, data = newdata)
                        frm <- formula(object$call)
			mt <- delete.response(terms(frm, data=newdata))
#			mf <- lm(object$formula, newdata, method="model.frame")
			mf <- model.frame(mt, newdata)
			X <- model.matrix(mt, mf)
			K <- ifelse(colnames(X)[1] == "(Intercept)", 2, 1)
			m <- ncol(X)
			WX <- matrix(nrow=nrow(X),ncol=(m-(K-1)))
			for (k in K:m) {
				wx <- lag.listw(listw, X[,k], 
				    zero.policy=zero.policy)
				if (any(is.na(wx))) 
				    stop("NAs in lagged independent variable")
				WX[,(k-(K-1))] <- wx
			}
			X <- cbind(X, WX)
#  accommodate aliased coefficients 120314
                        if (any(object$aliased))
                            X <- X[,-which(object$aliased)]
			trend <- X %*% B
                        if (power) {
                            W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
                            res <- c(as(powerWeights(W, rho=object$rho,
                                X=trend, order=order, tol=tol), "matrix"))
                        } else {
                            res <- c(invIrW(listw, object$rho) %*% trend)
                        }
                        if (legacy) {
			    signal <- object$rho * lag.listw(listw, 
				res, zero.policy=zero.policy)
			    res <- c(trend + signal)
                        } else {
                            signal <- res - trend
                        }
                        if (pred.se) {
                            samples <- attr(lagImpact, "samples")$samples
                            irho <- attr(lagImpact, "samples")$irho
                            drop2beta <- attr(lagImpact, "samples")$drop2beta
                            nSim <- nrow(samples)
                            outmat <- matrix(NA, ncol=nSim, nrow=nrow(X))
                            for (i in 1:nSim) {
                                B <- samples[i, -drop2beta]
                                trend <- X %*% B
                                rho <- samples[i, irho]
                                if (power) {
                                    res <- c(as(powerWeights(W, rho=rho,
                                    X=trend, order=order, tol=tol), "matrix"))
                                } else {
                                    res <- c(invIrW(listw, rho) %*% trend)
                                }
                                outmat[,i] <- res
                            }
                            pred.se <- apply(outmat, 1, sd)
                            attr(res, "pred.se") <- pred.se
                        }
			attr(res, "trend") <- c(trend)
			attr(res, "signal") <- c(signal)
		} else {
			if (is.null(listw) || !inherits(listw, "listw")) 
				stop ("spatial weights list required")
			if (nrow(newdata) != length(listw$neighbours))
				stop("mismatch between newdata and spatial weights")
			B <- object$coefficients
#			mt <- terms(object$formula, data = newdata)
                        frm <- formula(object$call)
			mt <- delete.response(terms(frm, data=newdata))
#			mt <- delete.response(terms(object$formula))
#			mf <- lm(object$formula, newdata, method="model.frame")
# resolved problem of missing response column in newdata reported by
# Christine N. Meynard, 060201
			mf <- model.frame(mt, newdata)
			if (dim(mf)[1] != length(listw$neighbours))
			    stop("missing values in newdata")
			X <- model.matrix(mt, mf)
#  accommodate aliased coefficients 120314
                        if (any(object$aliased))
                            X <- X[,-which(object$aliased)]
			trend <- X %*% B
                        if (power) {
                            W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
                            res <- c(as(powerWeights(W, rho=object$rho,
                                X=trend, order=order, tol=tol), "matrix"))
                        } else {
                            res <- c(invIrW(listw, object$rho) %*% trend)
                        }
                        if (legacy) {
			    signal <- object$rho * lag.listw(listw, 
				res, zero.policy=zero.policy)
			    res <- c(trend + signal)
                        } else {
                            signal <- res - trend
                        }
                        if (pred.se) {
                            samples <- attr(lagImpact, "samples")$samples
                            irho <- attr(lagImpact, "samples")$irho
                            drop2beta <- attr(lagImpact, "samples")$drop2beta
                            nSim <- nrow(samples)
                            outmat <- matrix(NA, ncol=nSim, nrow=nrow(X))
                            for (i in 1:nSim) {
                                B <- samples[i, -drop2beta]
                                trend <- X %*% B
                                rho <- samples[i, irho]
                                if (power) {
                                    res <- c(as(powerWeights(W, rho=rho,
                                    X=trend, order=order, tol=tol), "matrix"))
                                } else {
                                    res <- c(invIrW(listw, rho) %*% trend)
                                }
                                outmat[,i] <- res
                            }
                            pred.se <- apply(outmat, 1, sd)
                            attr(res, "pred.se") <- pred.se
                        }
			attr(res, "trend") <- c(trend)
			attr(res, "signal") <- c(signal)
		}
	}
	class(res) <- "sarlm.pred"
	res
}

print.sarlm.pred <- function(x, ...) {
	res <- data.frame(fit=as.vector(x), trend=attr(x, "trend"),
		signal=attr(x, "signal"))
	print(res, ...)
	invisible(res)
}

