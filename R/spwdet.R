# Copyright 2000-2 by Roger S. Bivand. 
#

spwdet <- function(sparseweights, rho)
{
	if(!inherits(sparseweights, "spatial.neighbour"))
             stop("Not a sparse weights object")
	if(missing(rho) || !is.numeric(rho))
		stop("rho incorrectly specified")
	n <- length(attr(sparseweights, "region.id"))
	size <- length(sparseweights$weights)
	vals <- -rho*sparseweights$weights
	z <- .C("spRdet",
			n = as.integer(n),
			size = as.integer(size),
			p1 = as.integer(sparseweights$from),
			p2 = as.integer(sparseweights$to),
			value = as.double(vals),
			determinant = double(1),
			pideterminant = double(1),
			exponent = integer(1),
			PACKAGE="spdep"
	)
	list(det=z$determinant, exp=z$exponent)
}

log.spwdet <- function(sparseweights, rho)
{
	if(!inherits(sparseweights, "spatial.neighbour"))
             stop("Not a sparse weights object")
	if(missing(rho) || !is.numeric(rho) ||
		rho >= 1 || rho <= -1)
		stop("rho incorrectly specified")
	res <- spwdet(sparseweights, rho=rho)
	det <- res$det * 10^res$exp
	if (det < .Machine$double.eps^0.5) {
		warning(paste("failure in spwdet:", det, res$det, res$exp))
		fres <- NaN
	} else	fres <- log(det)
	fres
}

