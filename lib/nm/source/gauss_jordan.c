#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
/////////////////////////////////////////////////////////////////////////
/// FIX-ME
/// Andrew, the following is going to cause trouble across elmer/elmer2.
/// If you need to know about cells, this should be in with mb_cns/mbcns2.
#include "../../../app/mb_cns/source/cns_cell.h"
/////////////////////////////////////////////////////////////////////////
#include "./gauss_jordan.h"
#define NR_END 1
#define FREE_ARG char*
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}


void nrerror(char error_text[]){
/* Numerical Recipes standard error handler */
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

/* void free_ivector(int *v, long nl, long nh) */
/* /\* free an int vector allocated with ivector() *\/ */
/* { */
/* 	free((FREE_ARG) (v+nl-NR_END)); */
/* } */

void free_ivector(int *v)
/* free an int vector allocated with ivector() */
{
	free(v);
}

void gaussj(struct cell_center *cell, int n, int m)
     /* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
      * is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
      * output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
      * vectors. 
      */
{    
    /*                    i  
     *                  ---->
     *         [                     ]
     *         [                     ]  
     *      |  [                     ]
     *   j  |  [                     ] 
     *      v  [                     ]
     *         [                     ]
     *         [                     ]
     *
     */

    int pivot_index_i, pivot_index_j;
    int i,j, row, col;
    float dum,pivinv;

    for (j = 1; j <= 5; j++) {
	for (i = 1; i <= 5; i++ ) {
	    cell->B_L[i][j] = cell->M_L[i][j];
	}
    }

    for (j = 1; j <=5; j++) {
	cell->h_L[1][j] = cell->r_L[1][j];
    }

    for (j=1;j<=n;j++) { /* This is the main loop over the rows.*/

	pivot_index_i = j;
	pivot_index_j = j;

	/* We are now ready to divide the pivot row by the */
	/* pivot element, located at irow and icol.*/

	if (cell->B_L[pivot_index_i][pivot_index_j] == 0.0) nrerror("gaussj: Singular Matrix");

	pivinv=1.0/cell->B_L[pivot_index_i][pivot_index_j];


	for (i=1;i<=n;i++) {
	    cell->B_L[i][pivot_index_j] *= pivinv;
	}

	for (i=1;i<=5;i++) {
	    cell->h_L[pivot_index_j][i] *= pivinv;
	}

	for (row=1;row<=n;row++) {/*  Next, we reduce the rows...*/
	    if (row != pivot_index_j) { /* ...except for the pivot one, of course.*/
		dum=cell->B_L[pivot_index_i][row];

		for (col=1;col<=n;col++) {
		    cell->B_L[col][row] -= cell->B_L[col][pivot_index_j]*dum;
		}
		for (col=1;col<=1;col++) {
		    cell->h_L[col][row] -= cell->h_L[col][pivot_index_j]*dum;
		}

	    }
	}
    }
    for (j = 1; j <=5; j++) {
	cell->r_L[1][j] = cell->h_L[1][j];
    }

}

void matrix_inverse(struct cell_center *cell, int n)
     /* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..2 x n][1..n]
      * is the input matrix with identity matrix on RHS. On
      * output, the RHS of a is the inverse
      */
{    
    /*                    i  
     *                  ---->
     *         [                     1    0    0    .....]
     *         [                     0    1    0    .....]  
     *      |  [                     0    0    1    .....]
     *   j  |  [                     0    0    0    .....] 
     *      v  [                     0    0    0    .....]
     *         [                     0    0    0    .....]
     *         [                     0    0    0    .....]
     *
     */

    int column, row, sub_row;
    double inverse, factor;

    for (row = 1; row <= n; ++row) {
	/* begin reducing the diagonal element and corresponding row */
	if(cell->M_L[row][row] == 0) {
	    printf("shit.... divide by zero\n");
	}
	inverse = 1/cell->M_L[row][row];

	for (column = 1; column <= 2*n; ++column) {
	    cell->M_L[column][row] *= inverse;
	}

	/* now the diagonal is 1 */
	
	for (sub_row = 1; sub_row <= n; ++sub_row) {
	    if (sub_row != row) {
		if (cell->M_L[row][sub_row] != 0) {
		    factor = cell->M_L[row][sub_row];
		    
		    for (column = 1; column <= 2*n; ++column) {
			cell->M_L[column][sub_row] = cell->M_L[column][sub_row] - factor * cell->M_L[column][row];
		    }
		    
		}
	    }
	}
	
    }
    
}

void inverse_product(struct cell_center *cell)
/* Takes inverse stored in M_L and multiplies by r_l */
{
    int j;
    double solution[5];

    for (j = 1; j <= 5; ++j) {
	solution[j] = cell->M_L[0][j] * cell->r_L[1][0] + cell->M_L[1][j] * cell->r_L[1][1] + cell->M_L[2][j] * cell->r_L[1][2] + 
	    cell->M_L[3][j] * cell->r_L[1][3] + cell->M_L[4][j] * cell->r_L[1][4];
    }

    for (j = 0; j < 5; ++j) {
	cell->r_L[1][j] = solution[j];
    }

}

void gaussj_pivoting(double **a, int n, double **b, int m)
     /* Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
      * is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
      * output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
      * vectors. 
      */
{    
    /*                    i  
     *                  ---->
     *         [                     ]
     *         [                     ]  
     *      |  [                     ]
     *   j  |  [                     ] 
     *      v  [                     ]
     *         [                     ]
     *         [                     ]
     *
     */

    int *indxc,*indxr,*ipiv;
    int i,icol,irow,j,k,l, row, col;
    float big,dum,pivinv,temp;
    indxc=ivector(1,n); /* The integer arrays ipiv, indxr, and indxc are 
			 * used for bookkeeping on the pivotiindxr=ivector(1,n); ng. */
    indxr=ivector(1,n);
    ipiv=ivector(1,n);
    for (j=1;j<=n;j++) ipiv[j]=0;
    for (i=1;i<=n;i++) { /* This is the main loop over the columns to be big=0.0; reduced.*/
	big = 0;
	for (j=1;j<=n;j++) {/*This is the outer loop of the search for a pivot if (ipiv[j] != 1) element.*/
	    if (ipiv[i] == 0) {
		if (fabs(a[i][j]) >= big) {
		    big=fabs(a[i][j]);
		    irow=j;
		    icol=i;
		}
	    }
	}
	printf("biggest element = %e at %i %i\n", big, irow, icol);
	++(ipiv[icol]);
	/* We now have the pivot element, so we interchange rows, if needed, to put the pivot
	 * element on the diagonal. The columns are not physically interchanged, only relabeled:
	 * indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
	 * indxr[i] is the row in which that pivot element was originally located. If indxr[i] =
	 * indxc[i] there is an implied column interchange. With this form of bookkeeping, the
	 * solution b will end up in the correct order, and the inverse matrix will be scrambled
	 * by columns.*/

	/* If maximum isn't on the diagonal, swap rows to make this so */
	if (irow != icol) {
	    for (l=1;l<=n;l++) SWAP(a[l][irow],a[l][icol])
	    for (l=1;l<=m;l++) SWAP(b[l][irow],b[l][icol])
	}

	indxr[i]=irow; /* We are now ready to divide the pivot row by the */
	indxc[i]=icol; /* pivot element, located at irow and icol.*/

	if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
	pivinv=1.0/a[icol][icol];

	for (l=1;l<=n;l++) {
	    a[l][icol] *= pivinv;
	    printf("%g ", a[l][icol]);
	}
       
	for (l=1;l<=m;l++) {
	    b[l][icol] *= pivinv;
	    printf("| %g\n", b[l][icol]);
	}

	for (row=1;row<=n;row++) {/*  Next, we reduce the rows...*/
	    if (row != icol) { /* ...except for the pivot one, of course.*/
		dum=a[icol][row];

		for (col=1;col<=n;col++) {
		    a[col][row] -= a[col][icol]*dum;
		    printf("%g ", a[col][row]);
		}
		for (col=1;col<=m;col++) {
		    b[col][row] -= b[col][icol]*dum;
		    printf("| %g\n", b[col][row]);
		}

	    }
	}
    }
    /* This is the end of the main loop over columns of the reduction. It only remains to unscramble
     * the solution in view of the column interchanges. We do this by interchanging pairs of
     * columns in the reverse order that the permutation was built up. */
    for (l=n;l>=1;l--) {
	if (indxr[l] != indxc[l])
	    for (k=1;k<=n;k++)
		SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    } /* And we are done.*/
    free_ivector(ipiv);
    free_ivector(indxr);
    free_ivector(indxc);
}
