# Copyright 2001-7 by Roger Bivand 
#

.spChkOption <- new.env(FALSE, globalenv())
assign("spChkID", FALSE, env = .spChkOption)
.conflicts.OK <- TRUE

.onLoad <- function(lib, pkg) {
	require(methods)
}

#.onLoad <- function(pkg, lib) {
#cat("spdep: a package for analysing spatial dependence\n")
#require(maptools)
#.First.lib <- function(lib, pkg) {
#	library.dynam("spdep", pkg, lib)
#}
#.noGenerics <- TRUE

