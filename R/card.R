# Copyright 2001 by Roger Bivand 
#

card <- function(nb) {
    if (class(nb) != "nb") stop("not a neighbours list")
    z <- .Call("card", nb, PACKAGE="spdep")
    invisible(z)
}
