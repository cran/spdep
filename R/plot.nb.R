# Copyright 2001 by Roger Bivand
#


plot.nb <- function(x, coords, col="black", points=TRUE, add=FALSE, 
	arrows=FALSE, length=0.1, ...) {
	nb <- x
	sym <- is.symmetric.nb(nb, verbose = FALSE, force = FALSE)
	x <- coords[,1]
	y <- coords[,2]
	n <- length(nb)
	xlim <- range(x)
	ylim <- range(y)
	if (!add) {
		plot.new()
        	plot.window(xlim = xlim, ylim = ylim, log="", asp=1)
	}
	cardnb <- card(nb)
	if (length(col) < n) col <- rep(col[1], n)
	for (i in 1:n) {
		if (cardnb[i] > 0) {
        		inb <- nb[[i]]
        		for (j in inb) {
				if (sym) {
					lines(c(x[i], x[j]), c(y[i], y[j]),
						col=col[i], ...)
				} else {
					if (arrows) 
						arrows(x[i], y[i], x[j], y[j], 
						col=col[i], length=length, ...)
					else lines(c(x[i], x[j]), c(y[i], y[j]),
						col=col[i], ...)
				}

			}
		}
	}
	if (points) points(x, y, ...)
}
