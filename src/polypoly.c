/* Copyright 2004 by Roger S. Bivand. */

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

SEXP polypoly(SEXP p1, SEXP n01, SEXP p2, SEXP n02, SEXP snap)
{
	int n1=INTEGER_POINTER(n01)[0], n2=INTEGER_POINTER(n02)[0], pc=0;
	int i, j, k=0;
	double sn=NUMERIC_POINTER(snap)[0], dist;
	double x1, x2, y1, y2;

	SEXP ans;
	PROTECT(ans = NEW_INTEGER(1)); pc++;

	for (i=0; (i < n1) && (k < 2); i++) {
		x1 = NUMERIC_POINTER(p1)[i];
		y1 = NUMERIC_POINTER(p1)[n1 + i];
		for (j=0; (j < n2) && (k < 2); j++) {
			x2 = NUMERIC_POINTER(p2)[j];
			y2 = NUMERIC_POINTER(p2)[n2 + j];
			dist = pythag((x1-x2), (y1-y2));
			if (dist < sn) k++;
			if (k > 1) break;
		}
	}
	
	INTEGER_POINTER(ans)[0] = k;

	UNPROTECT(pc); /* ans */
	return(ans);
}

