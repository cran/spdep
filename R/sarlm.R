# Copyright 2002 by Roger Bivand
#

residuals.sarlm <- function(object, ...) {
	object$residuals
}

deviance.sarlm <- function(object, ...) {
	deviance(object$lm.target)
}

coef.sarlm <- function(object, ...) {
	ret <- object$coefficients
	if(object$type == "error") ret <- c(ret, object$lambda)
	else ret <- c(ret, object$rho)
	ret
}

fitted.values.sarlm <- function(object, ...) {
	object$fitted.values
}

predict.sarlm <- function(object, newdata=NULL, listw=NULL, 
	zero.policy=FALSE, ...) {
	if (is.null(newdata)) {
		res <- fitted.values(object)
		X <- model.matrix(terms(object$lm.model), 
				model.frame(object$lm.model))
		B <- coefficients(object$lm.target)
		y <- model.response(model.frame(object$lm.model))
		tarX <- model.matrix(terms(object$lm.target), 
				model.frame(object$lm.target))
		tary <- model.response(model.frame(object$lm.target))
		if (object$type == "error") {
			attr(res, "trend") <- as.vector(X %*% B)
			attr(res, "signal") <- as.vector( -1 * (tary - y) - 					-1 * (tarX - X) %*% B)
		} else {
			attr(res, "trend") <- fitted.values(object$lm.target)
			attr(res, "signal") <- as.vector( -1 * (tary - y))
		}
	}
	else {
		if (is.null(listw) || !inherits(listw, "listw")) 
			stop ("spatial weights list required")
		if (nrow(newdata) != length(listw$neighbours))
			stop("mismatch between newdata and spatial weights")
		if (object$type == "error") {
			B <- coefficients(object$lm.target)
			tt <- terms(object$lm.model) 
			X <- model.matrix(delete.response(tt), data=newdata)
			trend <- X %*% B
			signal <- rep(0, length(trend))
			res <- trend + signal
			attr(res, "trend") <- trend
			attr(res, "signal") <- signal
		} else if (object$type == "mixed") {
			B <- coefficients(object$lm.target)
			mt <- terms(object$formula, data = newdata)
			mf <- lm(object$formula, newdata, method="model.frame")
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
			trend <- X %*% B
			signal <- object$rho * lag.listw(listw, 
				(invIrW(listw, object$rho) %*% trend), 
				zero.policy=zero.policy)
			res <- trend + signal
			attr(res, "trend") <- trend
			attr(res, "signal") <- signal

		} else {
			B <- coefficients(object$lm.target)
			mt <- terms(object$formula, data = newdata)
			mf <- lm(object$formula, newdata, method="model.frame")
			X <- model.matrix(mt, mf)
			trend <- X %*% B
			raw.sig <- invIrW(listw, object$rho) %*% trend
			signal <- object$rho * lag.listw(listw, 
				raw.sig, zero.policy=zero.policy)
			res <- trend + signal
			attr(res, "trend") <- trend
			attr(res, "signal") <- signal
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

