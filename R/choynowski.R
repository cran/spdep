# Copyright 2004 by Roger Bivand 
#

choynowski <- function(n, x, row.names=NULL, tol=.Machine$double.eps^0.5) {
  len <- length(n)
  res <- numeric(len)
  nsum <- sum(n)
  xsum <- sum(x)
  b <- nsum/xsum
  E <- x*b
  type <- (n < E)
  for (i in 1:len) {
    if(type[i]) {
      for (j in 0:n[i]) {
        xx <- (E[i]^j*exp(-E[i])) / gamma(j + 1)
        res[i] <- res[i] + xx
      }
    } else {
      xx <- 1
      x <- n[i]
      while (xx > tol) {
        xx <- (E[i]^x*exp(-E[i])) / gamma(x + 1)
        res[i] <- res[i] + xx
        x <- x + 1
      }
    }
  }
  if (is.null(row.names)) 
    res <- data.frame(pmap=res, type=type)
  else
    res <- data.frame(pmap=res, type=type, row.names=row.names)
  res
}
