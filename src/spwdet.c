/* Copyright 2000-2 by Roger S. Bivand. */

#include <stdio.h>
#include <math.h>
 #include <R.h>
#include <Rdefines.h>
#define ROFFSET 0 
#include "spMatrix.h"


void spRdet(int *n, int *size, int *p1, int *p2, double *value,
	double *determinant, double *piDeterminant, int *exponent)
{

	char *matrix;
	int i, errnum;
	const char *errormess[] = {
		"spOKAY",
		"spSMALL_PIVOT",
		"spZERO_DIAG",
		"spSINGULAR",
		"spNO_MEMORY",
       		"spPANIC"};
	
	matrix = spCreate(*n, 0, &errnum);
	if (errnum != 0) error ("error creating matrix: %s\n",
			errormess[errnum]);
	spClear(matrix);
	for (i=0; i<*size; i++) {
		spADD_REAL_ELEMENT(spGetElement(matrix,
			p1[i]-ROFFSET, p2[i]-ROFFSET), value[i]);
	}
	errnum = spError(matrix);
	if (errnum != 0) error ("error filling matrix with -rho*W: %s\n",
			errormess[errnum]);
	for (i=1; i<=*n; i++) {
		spADD_REAL_ELEMENT(spGetElement(matrix,i,i),
			1.0);
	}
	errnum = spError(matrix);
	if (errnum != 0) error ("error filling matrix with I: %s\n",
			errormess[errnum]);
	errnum = spOrderAndFactor(matrix, NULL, 0, 0, 1);
	if (errnum != 0) { if (errnum == 1)
		warning ("error ordering and factoring matrix: %s\n",
			errormess[errnum]);
		else {
			spDestroy(matrix);
			error ("error ordering and factoring matrix: %s\n",
			errormess[errnum]);
		}
	}
	spDeterminant(matrix, exponent, determinant);
	spDestroy(matrix);
	return;
}

