# Copyright 2001 by Roger Bivand
#


plot.nb <- function(x, coords, col="black", points=TRUE, add=FALSE, ...) {
	nb <- x
	x <- coords[,1]
	y <- coords[,2]
	n <- length(nb)
	xlim <- range(x)
	ylim <- range(y)
	if (!add) {
		plot.new()
        	plot.window(xlim = xlim, ylim = ylim, log="", asp=1)
	}
	for (i in 1:n) {
        	inb <- nb[[i]]
        	for (j in inb)
			lines(c(x[i], x[j]), c(y[i], y[j]),
				col=col, ...)
	}
	if (points) points(x, y, ...)
}
