/*
 *  based on code taken from:
 *  class/class.c by W. N. Venables and B. D. Ripley  Copyright (C) 1994-9
 *  and written by Roger Bivand (C) 2001
 */

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

#define MAX_TIES 1000

SEXP
dnearneigh(SEXP din1, SEXP din2, SEXP pnte, SEXP p, SEXP test)
{
    int   j, k, kn, npat, nte, pdim, pc=0;
    int   pos[MAX_TIES];
    double dist, tmp, dn, dn0;
    SEXP ans;
    SEXP class;
    SEXP nbtype;
    SEXP dists;
    
    dn0 = NUMERIC_POINTER(din1)[0];
    dn = NUMERIC_POINTER(din2)[0];
    nte = INTEGER_POINTER(pnte)[0];
    pdim = INTEGER_POINTER(p)[0];
    PROTECT(ans = NEW_LIST(1)); pc++;
    PROTECT(dists = NEW_NUMERIC(2)); pc++;
    NUMERIC_POINTER(dists)[0] = dn0;
    NUMERIC_POINTER(dists)[1] = dn;
    SET_VECTOR_ELT(ans, 0, NEW_LIST(nte));
    PROTECT(class = NEW_CHARACTER(1)); pc++;
    PROTECT(nbtype = NEW_CHARACTER(1)); pc++;
    SET_STRING_ELT(class, 0, COPY_TO_USER_STRING("nb"));
    SET_STRING_ELT(nbtype, 0, COPY_TO_USER_STRING("distance"));
    setAttrib(VECTOR_ELT(ans, 0), R_ClassSymbol, class);
    setAttrib(VECTOR_ELT(ans, 0), install("nbtype"), nbtype);
    setAttrib(VECTOR_ELT(ans, 0), install("distances"), dists);
    dn0 = dn0*dn0;
    dn = dn*dn;
    for (npat = 0; npat < nte; npat++) {
	kn = 0;
	for (j = 0; j < nte; j++) {
	    if (j == npat) continue;
	    dist = 0.0;
	    for (k = 0; k < pdim; k++) {
		tmp = NUMERIC_POINTER(test)[npat + k * nte]
			- NUMERIC_POINTER(test)[j + k * nte];
		dist += tmp * tmp;
	    }
	    if (dist > dn0 && dist <= dn) {
		pos[kn] = j;
		if (++kn == MAX_TIES - 1)
		    error("too many ties in dnearneigh");
	    }
	}
	if (kn < 1) {
	    SET_VECTOR_ELT(VECTOR_ELT(ans, 0), npat, NEW_INTEGER(1));
	    INTEGER_POINTER(VECTOR_ELT(VECTOR_ELT(ans, 0), npat))[0] = 0;
	} else {
	    SET_VECTOR_ELT(VECTOR_ELT(ans, 0), npat, NEW_INTEGER(kn));
	    for (k = 0; k < kn; k++) {
	    	INTEGER_POINTER(VECTOR_ELT(VECTOR_ELT(ans, 0), npat))[k]
		    = pos[k]+ROFFSET;
	    }
	}
    }
    UNPROTECT(pc); 
    return(ans);
}


