/* Copyright 2000-2 by Roger S. Bivand. */

#include <stdio.h>
#include <math.h>
 #include <R.h>
#include <Rdefines.h>
#define ROFFSET 0 
#include "spMatrix.h"


void spRdet(int *n, int *size, int *debug, int *p1, int *p2, double *value,
	double *determinant, double *piDeterminant, int *exponent)
{

	char *matrix;
	char *file="sparsestats";
	int i, errnum;
	int blockhits;
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
	if (*debug == 1) {
		errnum = spFileStats(matrix, file,
			       	"Matrix created and cleared");
		Rprintf("Matrix created and cleared\n");
	}
	blockhits = spHowManyHits(matrix);
	for (i=0; i<*size; i++) {
		spADD_REAL_ELEMENT(spGetElement(matrix,
			p1[i]-ROFFSET, p2[i]-ROFFSET), value[i]);
	}
	errnum = spError(matrix);
	if (errnum != 0) error ("error filling matrix with -rho*W: %s\n",
			errormess[errnum]);
	if (*debug == 1) {
		errnum = spFileStats(matrix, file, "filled with -rho*W");
		Rprintf("After filling with -rho*W\n");
	}
	if (blockhits < spHowManyHits(matrix))
		error("Suspicious allocation of memory in sparse functions\n         ...bailing out! Save workspace and quit R!\n");
	for (i=1; i<=*n; i++) {
		spADD_REAL_ELEMENT(spGetElement(matrix,i,i),
			1.0);
	}
	errnum = spError(matrix);
	if (errnum != 0) error ("error filling matrix with I: %s\n",
			errormess[errnum]);
	if (*debug == 1) {
		errnum = spFileStats(matrix, file, "filled diagonal");
		Rprintf("After filling with I\n");
	}
	if (blockhits < spHowManyHits(matrix))
		error("Suspicious allocation of memory in sparse functions\n         ...bailing out! Save workspace and quit R!\n");
	errnum = spOrderAndFactor(matrix, NULL, 0.01, 0, 1);
	if (errnum != 0) { if (errnum == 1)
		warning ("error ordering and factoring matrix: %s\n",
			errormess[errnum]);
		else {
			spDestroy(matrix);
			error ("error ordering and factoring matrix: %s\n",
			errormess[errnum]);
		}
	}
	if (*debug == 1) {
		Rprintf("After spOrderAndFactor\n");
		errnum = spFileStats(matrix, file,
				"after ordering and factoring");
	}
	if (blockhits < spHowManyHits(matrix))
		error("Suspicious allocation of memory in sparse functions\n         ...bailing out! Save workspace and quit R!\n");
	spDeterminant(matrix, exponent, determinant);
	if (*debug == 1) {
		Rprintf("After spDeterminant: %d %f\n",
			       	*exponent, *determinant);
	}
	spDestroy(matrix);
	return;
}

