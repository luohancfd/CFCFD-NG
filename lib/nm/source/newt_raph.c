/** \file newt_raph.c
 * \ingroup nm
 *  \brief A Newton-Raphson method to find the zeros for a set of equations
 *
 * This module provides a Newton-Raphson method to find the zeros for
 * a set of equations.  The equation and its Jacobian must be defined by the
 * user.
 *
 * \date 14-Apr-04
 * \author Rowan J Gollan
 * \version 1.1
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vecmat.h"
#include "newt_raph.h"

char nrid[] = "$Id$";

struct nrMem *newt_raph_init(int ndim,
			     int (*fun) (int, double *, double *, void *),
			     int (*Jfun) (int, double *, double **, void *),
			     void *fdata, void *Jdata, double tol)
{

    struct nrMem *nr_mem;

    nr_mem = (struct nrMem *) calloc(1, sizeof(struct nrMem));

    nr_mem->ndim = ndim;
    nr_mem->f = fun;
    nr_mem->Jac = Jfun;
    nr_mem->fdata = fdata;
    nr_mem->Jdata = Jdata;
    nr_mem->tol = tol;

    /* Assign working space */
    nr_mem->G = mem_alloc_vector(ndim);
    nr_mem->minusG = mem_alloc_vector(ndim);
    nr_mem->yold = mem_alloc_vector(ndim);
    nr_mem->ynew = mem_alloc_vector(ndim);
    nr_mem->dely = mem_alloc_vector(ndim);
    nr_mem->dGdy = mem_alloc_sqmat(ndim);

    return (struct nrMem *) nr_mem;

}

int newt_raph(struct nrMem *nr_mem, double *y0, double *y)
{

#   define MAX_COUNT 10

    int ndim, i, count;
    double *G, *minusG, *dely, *yold, *ynew;
    double **dGdy;
    double tol;

    ndim = nr_mem->ndim;
    tol = nr_mem->tol;
    
    G = nr_mem->G;
    minusG = nr_mem->minusG;
    dely = nr_mem->dely;
    yold = nr_mem->yold;
    ynew = nr_mem->ynew;
    dGdy = nr_mem->dGdy;

    /* Evaluate G(y(0)) based on initial guess */
    nr_mem->f(ndim, y0, G, nr_mem->fdata);
    set_minusG(ndim, G, minusG);
    nr_mem->Jac(ndim, y0, dGdy, nr_mem->Jdata);
    copy_vector(y0, yold, ndim);


    for (count = 0; count < MAX_COUNT; count++) {
        msys_solve(dGdy, dely, minusG, ndim);
        set_ynew(ndim, yold, dely, ynew);

        if (test_tol(ndim, dely, tol) == 1)
            break;
        else {
            copy_vector(ynew, yold, ndim);
            nr_mem->f(ndim, yold, G, nr_mem->fdata);
            set_minusG(ndim, G, minusG);
            nr_mem->Jac(ndim, yold, dGdy, nr_mem->Jdata);
        }

    }

    if (count >= MAX_COUNT) {
        printf("Newton-Raphson method did not converge (in newt_raph.c)\n");
        printf("Iterations count = %d, tolerance = %e \n", count, tol);
        print_vector(ndim, y0, "y0");
        print_vector(ndim, y, "y");

        return 1;
    }

    for (i = 0; i < ndim; i++) {
    	y[i] = ynew[i];
    }

#ifdef TEST_NEWT_RAPH
    printf("Iterations count in Newton-Raphson method: %d \n", count);
#endif

#undef MAX_COUNT

    return 0;
}

int set_minusG(int ndim, double *G, double *minusG)
{
    int i;
    for (i = 0; i < ndim; i++)
        minusG[i] = -1.0 * G[i];

    return 0;

}

int set_ynew(int ndim, double *yold, double *dely, double *ynew)
{
    int i;

    for (i = 0; i < ndim; i++)
        ynew[i] = yold[i] + dely[i];

    return 0;

}

int test_tol(int ndim, double *dely, double tol)
{

    int i, flag;

    flag = 1;   /* Set flag to pass test */
    for (i = 0; i < ndim; i++) {
	if (fabs(dely[i]) >= tol)
            flag = 0;   /* Set flag to fail if just one fails. */
    }

    return flag;
}

#ifdef TEST_NEWT_RAPH

int f1(int ndim, double *y, double *G, void *fdata)
{

    G[0] = 4.0 - y[0]*y[0] - y[1]*y[1];
    G[1] = 1.0 - exp(y[0]) - y[1];

    return 0;

}

int dfdy1(int ndim, double *y, double **dGdy, void *Jdata)
{

    dGdy[0][0] = -2.0 * y[0];   dGdy[0][1] = -2.0 * y[1];
    dGdy[1][0] = -1.0 * exp(y[0]); dGdy[1][1] = -1.0;

    return 0;

}


int f2(int ndim, double *y, double *G, void *fdata)
{

    G[0] = exp(y[0]) - y[1];
    G[1] = y[0] * y[1] - exp(y[0]);

    return 0;

}

int dfdy2(int ndim, double *y, double **dGdy, void *Jdata)
{

    dGdy[0][0] = exp(y[0]);   dGdy[0][1] = -1.0;
    dGdy[1][0] = y[1] - exp(y[0]); dGdy[1][1] = y[0];

    return 0;

}
int f3(int ndim, double *y, double *G, void *fdata)
{

    G[0] = y[0] * y[0] - 2.0 * y[0] - y[1] + 0.5;
    G[1] = y[0] * y[0] + 4.0 * y[1] * y[1] - 4.0;

    return 0;

}

int dfdy3(int ndim, double *y, double **dGdy, void *Jdata)
{

    dGdy[0][0] = 2.0 * y[0] - 2.0;   dGdy[0][1] = -1.0;
    dGdy[1][0] = 2.0 * y[0];         dGdy[1][1] = 8.0 * y[1];

    return 0;

}

int main()
{

    double *yguess, *y;
    double tol;
    int ndim;
    struct nrMem *nrP;

    ndim = 2;

    yguess = mem_alloc_vector(ndim);
    y = mem_alloc_vector(ndim);

    printf("======================================================\n");
    printf(" Test for: newt_raph.c                                \n");
    printf("======================================================\n");
    printf("%s\n\n", nrid);

    printf("--- Test case 1: 2 x 2 Nonlinear system ---\n");
    printf("\n");

    printf("This is from page 175 of Gerald and Wheatley.\n");
    printf("testing :: int newt_raph(int ndim, double *y0, double *y,\n");
    printf("                         int (*fun) (int, double *, double *, void *),\n");
    printf("                         int (*Jfun) (int, double *, double **, void *),\n");
    printf("                         void *fdata, void *Jdata, double tol)\n\n");

    printf("Solve for x and y the set of equations:\n\n");
    printf("x^2 + y^2 = 4 \n");
    printf("e^x + y = 1 \n");

    yguess[0] = 1.0; yguess[1] = -1.7;
    tol = 1.0e-6;

    printf("We take as a guess x = 1.0, y = -1.7 and specify tol = %g .\n", tol);
    print_vector(ndim, yguess, "guess");

    nrP = newt_raph_init(ndim, f1, dfdy1, NULL, NULL, tol);
    newt_raph(nrP, yguess, y);

    printf("After the call to newt_raph() :\n");
    print_vector(ndim, y, "y");

    printf("The answer given in Gerald and Wheatley is:\n");
    printf("x = 1.004167, y = -1.729635\n");
    printf("DONE\n\n");

    printf("--- Test case 2: 2 x 2 Nonlinear system ---\n");
    printf("\n");

    printf("This is from page 177 of Gerald and Wheatley.\n");

    printf("Solve for x and y the set of equations:\n\n");
    printf("e^x - y   = 0 \n");
    printf("x y - e^x = 0 \n");

    yguess[0] = 0.95; yguess[1] = 2.7;
    tol = 1.0e-7;

    printf("We take as a guess x = 0.95, y = 2.7 and specify tol = %g .\n", tol);
    print_vector(ndim, yguess, "guess");

    nrP = newt_raph_init(ndim, f2, dfdy2, NULL, NULL, tol);
    newt_raph(nrP, yguess, y);

    printf("After the call to newt_raph() :\n");
    print_vector(ndim, y, "y");

    printf("The answer given in Gerald and Wheatley is:\n");
    printf("x = 1.000000, y = 2.718282\n");
    printf("DONE\n\n");

    printf("--- Test case 3: 2 x 2 Nonlinear system ---\n");
    printf("\n");

    printf("This is Example 2.20 from Mathews.\n");

    printf("Solve for x and y the set of equations:\n\n");
    printf("x^2 - 2x - y + 0.5 \n");
    printf("x^2 + 4y^2 - 4 = 0 \n");

    yguess[0] = 2.00; yguess[1] = 0.25;
    tol = 1.0e-7;

    printf("We take as a guess x = 2.00, y = 0.25 and specify tol = %g .\n", tol);
    print_vector(ndim, yguess, "guess");

    nrP = newt_raph_init(ndim, f3, dfdy3, NULL, NULL, tol);
    newt_raph(nrP, yguess, y);

    printf("After the call to newt_raph() :\n");
    print_vector(ndim, y, "y");

    printf("The answer given in Mathews is:\n");
    printf("x = 1.900677, y = 0.311219\n");
    printf("DONE\n\n");

    printf("======================================================\n");
    printf(" End of test for: newt_raph.c                         \n");
    printf("======================================================\n");


    return 0;

}

#endif
