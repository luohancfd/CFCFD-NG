/** \file ode_be.c
 * \ingroup nm
 * \brief A backward Euler ODE solver.
 *
 * This module provides a numerical method to solve
 * an ODE in N-dimensions based on the backward Euler
 * method which is an implicit method.
 *
 * \date 14-Apr-04
 * \author Rowan J Gollan
 * \version 1.0
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "vecmat.h"
#include "newt_raph.h"
#include "ode_solvers.h"

char odbeid[] = "$Id";

#ifdef TEST_MAIN

double lambda1 = -0.5;
double lambda2 = -5.0;

#endif

int ode_be(int ndim, int (*ode_fun) (int , double *, double *, void *),
	   int (*ode_Jfun) (int, double *, double **, void *),
	   void *fdata, void *Jdata, double *yold, double *ynew, double h)
{

#   define TOL 1.0e-7

    struct ode_data od;
    struct ode_Jdata odJ;

    /* Get the pertinent data into ode_data structure */

    od.y_init = yold; /* Alright to point, we're not changing it */
    od.h = h;
    od.fdata = fdata;
    od.ode_fun = ode_fun;

    odJ.Jdata = Jdata;
    odJ.ode_Jfun = ode_Jfun;
    odJ.h = h;

    newt_raph(ndim, yold, ynew, newt_fun, newt_Jfun, &(od), &(odJ), TOL);

#   undef TOL

    return 0;

}

int mem_alloc_ode_data(struct ode_data *od, int ndim)
{

    od->y_init = mem_alloc_vector(ndim);
    return 0;

}

int newt_fun(int ndim, double *y, double *G, void *newt_data)
{

    int i;
    double h;
    double *y_init, *ydot;
    int (*f) (int, double *, double *, void *);
    struct ode_data *od;
    void *fdata;


    od = (struct ode_data *) newt_data;

    y_init = od->y_init;
    h = od->h;
    f = od->ode_fun;
    fdata = od->fdata;


    /* Evalute all ydots */
    ydot = mem_alloc_vector(ndim);
    f(ndim, y, ydot, fdata);

    for (i = 0; i < ndim; i++) {
	G[i] = y[i] - y_init[i] - h * ydot[i];
    }

    free(ydot);

    return 0;

}

int newt_Jfun(int ndim, double *y, double **dGdy, void *newt_Jdata)
{

    int (*Jfun) (int, double *, double **, void*);
    void *Jdata;
    double **hdfdy, **dfdy, **I;
    struct ode_Jdata *odJ;
    double h;


    hdfdy = mem_alloc_sqmat(ndim);
    dfdy = mem_alloc_sqmat(ndim);
    I = mem_alloc_sqmat(ndim);

    odJ = (struct ode_Jdata *) newt_Jdata;
    Jfun = odJ->ode_Jfun;
    Jdata = odJ->Jdata;
    h = odJ->h;

    Jfun(ndim, y, dfdy, Jdata);
    eye(I, ndim);
    msys_scalmult(h, dfdy, hdfdy, ndim);
    msys_subtract(dGdy, I, hdfdy, ndim);

    free_sqmat(hdfdy, ndim);
    free_sqmat(dfdy, ndim);
    free_sqmat(I, ndim);

    return 0;

}

#ifdef TEST_MAIN

int ode_f1(int ndim, double *y, double *ydot, void *ode_data)
{

    ydot[0] = lambda1 * y[0];
    ydot[1] = y[0] + lambda2 * y[1];

    return 0;

}

int Jode_f1(int ndim, double *y, double **dfdy, void *Jdata)
{

    dfdy[0][0] = lambda1;  dfdy[0][1] = 0.0;
    dfdy[1][0] = 1.0;      dfdy[1][1] = lambda2;

    return 0;

}

int ode_f2(int ndim, double *y, double *ydot, void *ode_data)
{

    ydot[0] = -0.04 * y[0] + 1.0e4 * y[1] * y[2];
    ydot[1] = 0.04 * y[0] - 1.0e4 * y[1] * y[2] - 3.0e7 * pow(y[1], 2);
    ydot[2] = 3.0e7 * pow(y[1], 2);

    return 0;

}

int Jode_f2(int ndim, double *y, double **dfdy, void *Jdata)
{

    dfdy[0][0] = -0.04;  dfdy[0][1] = 1.0e4 * y[2];  dfdy[0][2] = 1.0e4 * y[1];
    dfdy[1][0] = 0.04;   dfdy[1][1] = -1.0e4 * y[2] - 2.0 * 3.0e7 * y[1]; dfdy[1][2] = -1.0e4 * y[1];
    dfdy[2][0] = 0.0;    dfdy[2][1] = 2.0 * 3.0e7 * y[1]; dfdy[2][2] = 0.0;

    return 0;

}

int main()
{

    int ndim;
    double **A, *yold, *ynew, *anl;
    double dt, t, t_final, t_print, t_mult;

    printf("======================================================\n");
    printf(" Test for: ode_be.c                                \n");
    printf("======================================================\n");
    printf("%s\n\n", odbeid);

    printf("--- Test case 1: MECH3750 Assignment Question ---\n");
    printf("\n");

    ndim = 2;
    A = mem_alloc_sqmat(ndim);
    A[0][0] = lambda1;  A[0][1] = 0.0;
    A[1][0] = 1.0;      A[1][1] = lambda2;
    dt = 0.01;
    t = 0.0;
    t_final = 1.0;
    yold = mem_alloc_vector(ndim);
    ynew = mem_alloc_vector(ndim);
    yold[0] = 1.0; yold[1] = 0.0;

    printf("Solve the 2 x 2 system of ordinary differential equations:\n");
    printf("y' = A y where\n");
    print_sqmat(ndim, A, "A");
    printf("using timestep = %g and the backwards Euler method.\n\n", dt);

    printf("testing ::: int ode_be(int ndim,\n");
    printf("                       int (*ode_fun) (int , double *, double *, void *),\n");
    printf("                       int (*ode_Jfun) (int, double *, double **, void *),\n");
    printf("                       void *fdata, void *Jdata, double *yold, double *ynew,\n");
    printf("                       double h)\n\n");

    while ( t < t_final ) { 

        ode_be(ndim, ode_f1, Jode_f1, NULL, NULL, yold, ynew, dt);
	copy_vector(ynew, yold, ndim);
	t = t + dt;

    } 
    print_vector(ndim, yold, "y(x = 1.0)");
    printf("The analytical solution is..\n");
    anl = mem_alloc_vector(ndim);
    anl[0] = 0.6065306;  anl[1] = 0.1332873;
    print_vector(ndim, anl, "analytical");
    printf("DONE\n\n");
    free_sqmat(A, ndim);
    free(yold);
    free(ynew);


    printf("--- Test case 2: from CVODE example ---\n");
    printf("\n");

    ndim = 3;
    yold = mem_alloc_vector(ndim);
    ynew = mem_alloc_vector(ndim);
    t = 0.0;
    t_final = 4.0e1;
    dt = 0.0002;
    t_print = 4.0e-1;
    t_mult = 10.0;
    yold[0] = 1.0; yold[1] = 0.0; yold[3] = 0.0;

    printf("Solve the following chemical kinetics problem..\n");
    printf("dy1/dt = -0.04 * y1 + 1.0e4*y2*y3\n");
    printf("dy2/dt = 0.04 * y1 - 1.e4*y2*y3 - 3.0e7*(y2)^2\n");
    printf("dy3/dt = 3.0e7*(y2)^2\n\n");
    printf("on the interval t = 0.0 to t = 4.0e10 with\n\n");
    print_vector(ndim, yold, "y0");
    
    while ( t <= t_final) {
	ode_be(ndim, ode_f2, Jode_f2, NULL, NULL, yold, ynew, dt);
	copy_vector(ynew, yold, ndim);
	t = t + dt;
	
	if (fabs(t_print - t) < 1.0e-8) {
	    t_print *= t_mult;
	    printf("At t = %8.4e \t y =  %12.6e \t %12.6e \t %12.6e \n",
                    t, yold[0], yold[1], yold[2]);
	}

    }

    printf("\nThe solution given in the CVODE report is:\n\n");
    
    printf("At t = 4.0000e-01 \t y =  9.851641e-01 \t 3.386242e-05 \t 1.480205e-02 \n");
    printf("At t = 4.0000e+00 \t y =  9.055097e-01 \t 2.240338e-05 \t 9.446793e-02 \n");
    printf("At t = 4.0000e-01 \t y =  7.158016e-01 \t 9.185045e-06 \t 2.841893e-01 \n");

    printf("DONE\n\n");

    printf("======================================================\n");
    printf(" End of test for: ode_be.c                            \n");
    printf("======================================================\n");

	    

    return 0;

} /* end int main() */

#endif
