# Copyright 2001-10 by Roger Bivand 
#

.spdepOptions <- new.env(FALSE, globalenv())
assign("spChkID", FALSE, env = .spdepOptions)
assign("zeroPolicy", FALSE, env = .spdepOptions)
assign("verbose", FALSE, env = .spdepOptions)
assign("cl", NULL, env = .spdepOptions)
assign("rlecuyerSeed", rep(12345, 6), env = .spdepOptions)

#.conflicts.OK <- TRUE

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

