# Copyright 2001-2 by Roger Bivand 
#

.First.lib <- function(lib, pkg) {
	library.dynam("spdep", pkg, lib)
}
