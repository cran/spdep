/* Copyright 2001 by Roger S. Bivand. */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

SEXP gsymtest(SEXP nb, SEXP glist, SEXP card)
{
	int i, icard, j, k, k1, n=length(nb), pc=0;
	double g;
	SEXP ans;
	PROTECT(ans = NEW_LOGICAL(1)); pc++;
	LOGICAL_POINTER(ans)[0] = TRUE;

	for (i=0; i < n; i++) {
	    icard = INTEGER_POINTER(card)[i];
	    for (j=0; j<icard; j++) {
		k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j];
		g = NUMERIC_POINTER(VECTOR_ELT(glist, i))[j];
		if (k > 0 && k <= n) {
		    for (k1=0; k1<INTEGER_POINTER(card)[k-ROFFSET]; k1++) {
			if (i+ROFFSET == INTEGER_POINTER(VECTOR_ELT(nb,
			    k-ROFFSET))[k1]) {
			    if (g != NUMERIC_POINTER(VECTOR_ELT(glist,
			        k-ROFFSET))[k1]) {
				LOGICAL_POINTER(ans)[0] = FALSE;
				UNPROTECT(pc);
				return(ans);
			    }
			}
		    }
		}
	    }
	}

	UNPROTECT(pc); /* ans */
	return(ans);
}

