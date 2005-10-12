# Copyright 2001-5 by Roger Bivand 
#

.spChkOption <- new.env(FALSE, globalenv())
assign("spChkID", FALSE, env = .spChkOption)
#.onLoad <- function(pkg, lib) {
cat("spdep: a package for analysing spatial dependence\n")
#require(maptools)
#.First.lib <- function(lib, pkg) {
#	library.dynam("spdep", pkg, lib)
#}
.noGenerics <- TRUE

