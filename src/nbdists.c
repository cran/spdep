/* Copyright 2000 by Roger S. Bivand. */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

SEXP nbdists(SEXP nb, SEXP x, SEXP np, SEXP dim)
{
	int i, j, j1, k, m, n, d, pc=0;
	SEXP ans;
        SEXP class;
	double tmp, tmp1;
	
	PROTECT(ans = NEW_LIST(1)); pc++;
	n = INTEGER_POINTER(np)[0];
	SET_VECTOR_ELT(ans, 0, NEW_LIST(n));
	d = INTEGER_POINTER(dim)[0];
	PROTECT(class = NEW_CHARACTER(1)); pc++;
	SET_STRING_ELT(class, 0, COPY_TO_USER_STRING("nbdist"));
	setAttrib(VECTOR_ELT(ans, 0), R_ClassSymbol, class);

	for (i=0; i < n; i++) {
		k = length(VECTOR_ELT(nb, i));
		if (k == 1 && INTEGER_POINTER(VECTOR_ELT(nb, i))[0] == 0) {
			SET_VECTOR_ELT(VECTOR_ELT(ans, 0), i,
				NEW_NUMERIC(1));
			NUMERIC_POINTER(VECTOR_ELT(VECTOR_ELT(ans, 0), i))[0]
				= NA_REAL;
		} else {
			SET_VECTOR_ELT(VECTOR_ELT(ans, 0), i,
				NEW_NUMERIC(k));
			for (j=0; j < k; j++) {
				j1 = INTEGER_POINTER(VECTOR_ELT(nb, i))[j]
					- ROFFSET;
				tmp1 = 0;
				for (m=0; m < d; m++) {
					tmp = NUMERIC_POINTER(x)[i + m * n]
					- NUMERIC_POINTER(x)[j1 + m * n];
					tmp1 += tmp * tmp;
				}
				NUMERIC_POINTER(VECTOR_ELT(VECTOR_ELT(ans, 0),
					i))[j] = sqrt(tmp1);
			}
		}
	}
	UNPROTECT(pc);
	return(ans);
}

