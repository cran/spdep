# Copyright 2001-3 by Roger Bivand 
#

.spChkOption <- new.env(FALSE, globalenv())
assign("spChkID", FALSE, env = .spChkOption)
.DESC <- packageDescription("spdep")
.spdep.Version <- paste(.DESC[["Package"]], ", version ", .DESC[["Version"]],
	 ", ", .DESC[["Date"]], sep="")
.spdep.Build <- paste("build:", .DESC[["Built"]])
#.onLoad <- function(pkg, lib) {
cat(paste(spdep()[1], ":\n a package for analysing spatial dependence,\n",
	" use help(get.spChkOption) for help on integrity checking\n",
	sep=""))
#require(maptools)
#.First.lib <- function(lib, pkg) {
#	library.dynam("spdep", pkg, lib)
#}
.noGenerics <- TRUE

