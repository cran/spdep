# Copyright 2001-3 by Roger Bivand 
#

.spChkOption <- new.env(FALSE, globalenv())
assign("spChkID", FALSE, env = .spChkOption)
.DESC <- read.dcf(system.file("DESCRIPTION", package="spdep"))
.spdep.Version <- paste(.DESC[1,1], ", version ", .DESC[1,2], ", ", .DESC[1,3], 
	sep="")
.spdep.Build <- paste("build:", .DESC[1, dim(.DESC)[2]])
#.onLoad <- function(pkg, lib) {
cat(paste(spdep()[1], ":\n a package for analysing spatial dependence,\n",
	" use help(get.spChkOption) for help on integrity checking\n",
	sep=""))
#require(maptools)
#.First.lib <- function(lib, pkg) {
#	library.dynam("spdep", pkg, lib)
#}
.noGenerics <- TRUE

