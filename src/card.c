/* Copyright 2000-2 by Roger S. Bivand. */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

SEXP card(SEXP nb)
{
	int i, n=length(nb), pc=0;
	SEXP ans;
	PROTECT(ans = NEW_INTEGER(n)); pc++;

	for (i=0; i < n; i++) {
	    if (INTEGER_POINTER(VECTOR_ELT(nb, i))[0] == 0) 
		INTEGER_POINTER(ans)[i] = 0;
	    else
		INTEGER_POINTER(ans)[i] = length(VECTOR_ELT(nb, i));
	}

	UNPROTECT(pc); /* ans */
	return(ans);
}

