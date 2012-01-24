/** \file vecmat.c
 * \ingroup nm
 * \brief A module providing some vector and matix functions.
 *
 * This module conatins some basic service functions (memory allocation)
 * and mathematical routines required when working with vectors and matrices.
 *
 * \date 15-Apr-04
 * \author Rowan J Gollan
 * \version 1.0
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vecmat.h"

char svnid[] = "$Id$";


double *mem_alloc_vector(int ndim)
{

    double *G;

    G = (double *) calloc(ndim, sizeof(double));
    if (G == NULL) {
        printf("Could not allocate memory for vector in vecmat.c \n");
        exit(-1);
    }

    return G;

}

int *mem_alloc_iarray(int ndim)
{

    int *G;
    G = (int *) calloc(ndim, sizeof(int));
    if (G == NULL) {
        printf("Could not allocate memory for vector in vecmat.c \n");
        exit(-1);
    }

    return G;

}


double **mem_alloc_sqmat(int ndim)
{

    double **A;

    int i, flag;
    flag = 0;

    A = (double **) calloc(ndim, sizeof(double *));
    if (A == NULL)
        ++flag;

    for (i = 0; i < ndim; i++) {
        A[i] = (double *) calloc(ndim, sizeof(double));
        if (A[i] == NULL)
            ++flag;
    }

    if (flag != 0) {
        printf("Could not allocate memory for square matrix in vecmat.c, flag = %d.\n", flag);
        exit(-1);
    }

    return A;

}

int free_sqmat(double **A, int ndim)
{

    int i;

    for(i = 0; i < ndim; i++) {
	free(A[i]);
    }

    free(A);

    return 0;

}

int copy_vector(double *a, double *b, int ndim)
{

   int i;

   for (i = 0; i < ndim; i++) {
       b[i] = a[i];
   }

   return 0;

}


int print_vector(int ndim, double *vec, char *name)
{

    int i;

    for (i = 0; i < ndim; i++) {
        printf("%s[%d] = %7.7g \n", name, i, vec[i]);
    }

    printf("\n");

    return 0;

}

int scale_vector(int ndim, double factor, double *v_in, double *v_out)
{
    int i;
    for (i = 0; i < ndim; i++) {
	v_out[i] = factor * v_in[i];
    }
    return (0);
}

int add_vectors(int ndim, double *a, double *b, double *c)
{
    int i;
    for (i = 0; i < ndim; ++i) {
	c[i] = a[i] + b[i];
    }
    return (0);
}

int subtract_vectors(int ndim, double *a, double *b, double *c)
{
    int i;
    for (i = 0; i < ndim; ++i) {
	c[i] = a[i] - b[i];
    }
    return (0);
}

int add_scalar_to_vector(int ndim, double scalar, double *a, double *b)
{
    int i;
    for (i = 0; i < ndim; ++i) {
	b[i] = scalar + a[i];
    }
    return (0);
}

int element_divide(int ndim, double *a, double *b, double *c)
{
    int i;
    for (i = 0; i < ndim; ++i) {
	c[i] = a[i] / b[i];
    }
    return (0);
}

int abs_vector(int ndim, double *a, double *b)
{
    int i;
    for (i = 0; i < ndim; ++i) {
	b[i] = fabs(a[i]);
    }
    return (0);
}

double max_vector(int ndim, double *a)
{
    int i;
    double max;
    max = a[0];

    for (i = 1; i < ndim; ++i) {
	if (a[i] > max)
	    max = a[i];
    }
    return (max);
}


int print_sqmat(int ndim, double **sqmat, char *name)
{

    int i, j;
   
    printf("%s = \n", name);

    for(i = 0; i < ndim; i++) {
	for(j = 0; j < ndim; j++) {
	    printf("%7.7g  ", sqmat[i][j]);
	}
	printf("\n");
    }

    printf("\n");

    return 0;

}


int eye(double **I, int ndim)
{

    int i, j;
    
    for (i = 0; i < ndim; i++) {
	for (j = 0; j < ndim; j++) {
	    if ( i == j )
		I[i][j] = 1.0;
	    else
		I[i][j] = 0.0;
	}
    }

    return 0;

}

int msys_subtract(double **A, double **B, double **C, int ndim)
{

    int i,j;

    for (i = 0; i < ndim; i++) {
	for (j = 0; j < ndim; j++) {
	    A[i][j] = B[i][j] - C[i][j];
	}
    }

    return 0;

}

int msys_scalmult(double C, double **A, double **B, int ndim)
{

    int i, j;
    
    for (i = 0; i < ndim; i++) {
	for (j = 0; j < ndim; j++) {
	    B[i][j] = C * A[i][j];
	}
    }

    return 0;

}
	    

/**
 * \brief Solve the matrix system A x = b.
 *
 * This function solves a matrix system of the form
 * \f$ \mathbf{A} \vec{x} = \vec{b} \f$ using the
 * Gaussian elimination method with partial pivoting.
 *
 * \b Inputs
 * - matrix A
 * - vector x: unknowns
 * - vector b: RHS
 *
 * \b Outputs
 * - updates vector x with solution
 *
 * \b Revisions
 * - 15-Apr-04 --- first attempt.
 *
 * \return int
 *
 * \author Rowan J Gollan (r.gollan@uq.edu.au)
 * \date 15-Apr-04
 * \version 1.0
 **/

int msys_solve(double **A, double *x, double *b, int ndim)
{

    int i, j, k, pv_temp;
    double pvt, tmp_val;
    
    double *tmp_row;
    int *pivot;

    pivot = mem_alloc_iarray(ndim);

    /* Elimination and partial-pivoting */

    for (j = 0; j < (ndim - 1); j++) {
	/* Find the value at present diagonal */
	pvt = fabs(A[j][j]);
	pivot[j] = j;
	pv_temp = j;

	/* Test the rows underneath for a larger pivot */
	for (i = (j + 1); i < ndim; i++) {
	    if (fabs(A[i][j]) > pvt) {
		pvt = fabs(A[i][j]);
		pv_temp = i;
	    }

	}

	/* See if we need to switch */
	if (pivot[j] != pv_temp) {
	    tmp_row = A[j];
	    A[j] = A[pv_temp];
	    A[pv_temp] = tmp_row;
	    tmp_val = b[j];
	    b[j] = b[pv_temp];
	    b[pv_temp] = tmp_val;
	}

	/* store multipliers */

	for (i = j + 1; i < ndim; i++) {
	    A[i][j] = A[i][j] / A[j][j];
	}

	/* Create zeros below main diagonal */


	for (i = j + 1; i < ndim; i++) {
	    for (k = j + 1; k < ndim; k++) {
		A[i][k] = A[i][k] - A[i][j] * A[j][k];
	    }
	    b[i] = b[i] - A[i][j] * b[j];
	}

    }
    /* Back substitution */
     
    x[ndim - 1] = b[ndim - 1] / A[ndim - 1][ndim - 1];

    for (j = (ndim - 2); j >= 0; j--) {
	x[j] = b[j];
	for (k = (ndim - 1); k >= j + 1; k--) {
	    x[j] = x[j] - x[k] * A[j][k];
	}
	x[j] = x[j] / A[j][j];
    }

    free(pivot);

    return 0;

}
    

#ifdef TEST_MAIN

int main()
{

    double *a, *b, *D, *x;
    double **A, **B, **C, **I, **F;
    int *m, *n;
    int size1, size2;
    char vecD[] = "D";
    char matF[] = "F";
    char matI[] = "I";

    size1 = 3;
    size2 = 4;


    printf("======================================================\n");
    printf(" Test for: vecmat.c                                   \n");
    printf("======================================================\n");
    printf("%s\n\n", svnid);

    printf("------------------------------------------------------\n");
    printf("     Service Function Tests\n");
    printf("------------------------------------------------------\n\n");

    printf("testing ::: int mem_alloc_vector(double *G, int ndim)\n");
    a = mem_alloc_vector(size1);
    printf("passed! \n\n");

    printf("testing ::: int mem_alloc_iarray(int *G, int ndim)\n");
    m = mem_alloc_iarray( size1 );
    printf("passed! \n\n");

    printf("testing ::: int mem_alloc_sqmat(double **A, int ndim) \n");
    A = mem_alloc_sqmat(size1);
    printf("passed! \n\n");

    printf("testing ::: int print_vector(int ndim, double *vec, char *name) \n\n");
    printf("We expect to see a vector \'D\' with entries (-1.0, 2.5, 3.1) \n");
    D = mem_alloc_vector(size1);
    D[0] = -1.0;
    D[1] = 2.5;
    D[2] = 3.1;
    print_vector(size1, D, vecD);
    printf("DONE \n\n");
 
    printf("testing ::: int print_sqmat(int ndim, double **sqmat, char *name) \n\n");
    printf("We expect to see a matrix \'F\' which is 3 x 3 and has entries\n");
    printf("[ [2.1, 5.8, 9778.0], [43.4476, -89.1, 1.002], [ 1000.0, -1000.0, 17.0] ] \n");
    F = mem_alloc_sqmat(size1);
    F[0][0] = 2.1;      F[0][1] = 5.8;       F[0][2] = 9778.0;
    F[1][0] = 43.4476;  F[1][1] = -89.1;     F[1][2] = 1.002;
    F[2][0] = 1000.0;   F[2][1] = -1000.0;   F[2][2] = 17.0;
    print_sqmat(size1, F, matF);
    printf("DONE \n\n");

    printf("testing ::: int eye(double **I, int ndim) \n\n");
    printf("We expect to see a 5 x 5 identity matrix.\n");
    printf("First we print \'I\' immediately after memory allocation\n");
    I = mem_alloc_sqmat(size2);
    print_sqmat(size2, I, matI);
    printf("Now we make the call to eye()\n");
    eye( I, size2);
    print_sqmat(size2, I, matI);
    printf("DONE\n\n");


    printf("------------------------------------------------------\n");
    printf("     Work-doing Function Tests\n");
    printf("------------------------------------------------------\n\n");

    printf("--- Test case 1: Matrix subtraction ---\n");
    printf("\n");

    printf("testing ::: int msys_subtract(double **A, double **B,\n");
    printf("                              double **C, int ndim)\n\n");
    printf("Find A = B - C where \n");
    B = mem_alloc_sqmat(size1);
    B[0][0] = 7.3; B[0][1] = 3.9; B[0][2] = 15.0;
    B[1][0] = 4.0; B[1][1] = -4.0; B[1][2] = -46.3;
    B[2][0] = 0.0; B[2][1] = 18.7; B[2][2] = 3.1;
    print_sqmat(size1, B, "B");
    C = mem_alloc_sqmat(size1);
    C[0][0] = 7.3; C[0][1] = -3.9; C[0][2] = 21.4;
    C[1][0] = 5.0; C[1][1] = -4.2; C[1][2] = -46.301;
    C[2][0] = 3.7; C[2][1] = 17.2; C[2][2] = 8.0;
    print_sqmat(size1, C, "C");
    printf("Therefore...\n");
    msys_subtract(A, B, C, size1);
    print_sqmat(size1, A, "A");
    printf("DONE\n\n");

    printf("--- Test case 2: Matrix system solve ---\n");
    printf("\n");
    printf("This is Example 2.2.2 in Schilling & Harris\n");
    printf("testing ::: int msys_solve(double **A, double *x,\n");
    printf("                           double *b, int ndim)\n\n");

    printf("Solve the system A x = b where\n");
    A[0][0] = 1.0; A[0][1] = -1.0; A[0][2] = 0.0;
    A[1][0] = -2.0; A[1][1] = 2.0; A[1][2] = -1.0;
    A[2][0] = 0.0; A[2][1] = 1.0; A[2][2] = -2.0;
    print_sqmat(size1, A, "A");
    printf("and\n");
    b = mem_alloc_vector(size1);
    b[0] = 2.0; b[1] = 1.0; b[2] = 6.0;
    print_vector(size1, b, "b");
    x = mem_alloc_vector(size1);

    msys_solve(A, x, b, size1);

    printf("After call to msys_solve() :\n");
    print_sqmat(size1, A, "A");
    print_vector(size1, b, "b");
    print_vector(size1, x, "x");

    printf("The answer should be x = [ 2.0, 0.0, -3.0 ] \n");
    printf("DONE\n\n");

    free(A); free(b); free(x);

    printf("--- Test case 3: Matrix system solve ---\n");
    printf("\n");
    printf("This is from page 130 of Gerald and Wheatley.\n");
    printf("testing ::: int msys_solve(double **A, double *x,\n");
    printf("                           double *b, int ndim)\n\n");

    printf("Solve the system A x = b where\n");
    A = mem_alloc_sqmat(size2);
    A[0][0] = 0.0; A[0][1] = 2.0; A[0][2] = 0.0; A[0][3] = 1.0;
    A[1][0] = 2.0; A[1][1] = 2.0; A[1][2] = 3.0; A[1][3] = 2.0;
    A[2][0] = 4.0; A[2][1] = -3.0; A[2][2] = 0.0; A[2][3] = 1.0;
    A[3][0] = 6.0; A[3][1] = 1.0; A[3][2] = -6.0; A[3][3] = -5.0;
    print_sqmat(size2, A, "A");
    printf("and\n");
    b = mem_alloc_vector(size2);
    b[0] = 0.0; b[1] = -2.0; b[2] = -7.0; b[3] = 6.0;
    print_vector(size2, b, "b");
    x = mem_alloc_vector(size2);

    msys_solve(A, x, b, size2);

    printf("After call to msys_solve() :\n");
    print_sqmat(size2, A, "A");
    print_vector(size2, b, "b");
    print_vector(size2, x, "x");

    printf("The answer should be x = [ -0.5, 1, 0.333, -2 ] \n");
    printf("DONE\n\n");

    printf("======================================================\n");
    printf(" End of test for: vecmat.c                            \n");
    printf("======================================================\n");

    return 0;

}

#endif

    
