# Copyright 2001-3 by Roger Bivand 
#

.spChkOption <- new.env(FALSE, globalenv())
assign("spChkID", FALSE, env = .spChkOption)
cat("Spatial object identity integrity check: FALSE\n")

.First.lib <- function(lib, pkg) {
	library.dynam("spdep", pkg, lib)
}
