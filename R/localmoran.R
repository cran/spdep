# Copyright 2001 by Roger Bivand 
#

localmoran <- function(x, listw, zero.policy=FALSE, spChk=NULL)
{
	if (class(listw) != "listw")
		stop(paste(deparse(substitute(listw)), "is not a listw object"))
	if (!is.null(attr(listw$neighbours, "self.included")) &&
		attr(listw$neighbours, "self.included"))
		stop("Self included among neighbours")
	n <- length(listw$neighbours)
	if (!is.numeric(x))
		stop(paste(deparse(substitute(x)), "is not a numeric vector"))
	if (any(is.na(x))) stop(paste("NA in ", deparse(substitute(x))))
	if (n != length(x))stop("Different numbers of observations")
	if (is.null(spChk)) spChk <- get.spChkOption()
	if (spChk && !chkIDs(x, listw))
		stop("Check of data and weights ID integrity failed")
	res <- data.frame(matrix(nrow=n,ncol=4))
	colnames(res) <- c("Ii", "E.Ii", "Var.Ii", "Z.Ii")
	z <- scale(x, scale=FALSE)
	lz <- lag.listw(listw, z, zero.policy=zero.policy)
	s2 <- sum(z^2)/n
	res[,1] <- (z/s2) * lz
	Wi <- sapply(listw$weights, sum)
	res[,2] <- -Wi / (n-1)
	b2 <- (sum(z^4)/n)/(s2^2)
	A <- (n-b2) / (n-1)
	B <- (2*b2 - n) / ((n-1)*(n-2))
	C <- Wi^2 / ((n-1)^2)
	Wi2 <- sapply(listw$weights, function(x) sum(x^2))
	Wikh2 <- sapply(listw$weights, function(x) {
		ifelse(is.null(x), 0, 1 - crossprod(x,x))
	})
	res[,3] <- A*Wi2 + B*Wikh2 - C
	res[,4] <- (res[,1] - res[,2]) / sqrt(res[,3])
	attr(res, "call") <- match.call()
	class(res) <- c(class(res), "localmoran")
	invisible(res)
}


