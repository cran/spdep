/* Copyright 2003 by Roger S. Bivand. */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

SEXP jcintern(SEXP nb, SEXP weights, SEXP dum, SEXP card) {
	int i, j, k, n=length(card), pc=0;
	double sum, sum1, wt;
	SEXP ans;
	PROTECT(ans = NEW_NUMERIC(1)); pc++;

	sum1 = 0.0;
	for (i=0; i < n; i++) {
	    sum = 0.0;
	    if (INTEGER_POINTER(card)[i] > 0) {
		for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
		    k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j];
		    wt = NUMERIC_POINTER(VECTOR_ELT(weights, i))[j];
		    sum += ((double) INTEGER_POINTER(dum)[k-ROFFSET]) * wt;
		}
		sum1 += ((double) INTEGER_POINTER(dum)[i]) * sum;
	    }
        }
	NUMERIC_POINTER(ans)[0] = sum1;
	
	UNPROTECT(pc); /* ans */
	return(ans);
}


/* SEXP jcdiffs(SEXP nb, SEXP weights, SEXP ifx, SEXP card) {
	int i, j, k, n=length(card), pc=0, iifx, kifx;
	double sum, sum1, wt;
	SEXP ans;
	PROTECT(ans = NEW_NUMERIC(1)); pc++;

	sum1 = 0.0;
	for (i=0; i < n; i++) {
	    sum = 0.0;
	    if (INTEGER_POINTER(card)[i] > 0) {
	    	iifx = INTEGER_POINTER(ifx)[i];
		for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
		    k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j];
		    kifx = INTEGER_POINTER(ifx)[k-ROFFSET];
		    wt = (double) NUMERIC_POINTER(VECTOR_ELT(weights, i))[j];
		    sum += (iifx == kifx) ? (double) 0 : wt;
		}
		sum1 += sum;
	    }
        }
	NUMERIC_POINTER(ans)[0] = sum1/2;
	
	UNPROTECT(pc); 
	return(ans);
}
*/

