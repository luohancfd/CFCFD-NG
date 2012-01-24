/** \file newt_raph.h
 * \ingroup nm
 * \brief Header file for newt_raph.c
 *
 * \date 25-Apr-04
 * \author Rowan J Gollan
 * \version 1.0
 *
 */

#ifndef NR_HEADER_ALREADY_INCLUDED

struct nrMem {

    int ndim;

    int (*f)   (int, double *, double *, void *);
    int (*Jac) (int, double *, double **, void *);

    void *fdata;
    void *Jdata;
    double tol;

    /* Working arrays */
    double *G;
    double *minusG;
    double *yold;
    double *ynew;
    double *dely;
    double **dGdy;

};


struct nrMem *newt_raph_init(int ndim,
			     int (*fun) (int, double *, double *, void *),
			     int (*Jfun) (int, double *, double **, void *),
			     void *fdata, void *Jdata, double tol);
int newt_raph(struct nrMem *nr_mem, double *y0, double *y);
int set_minusG(int ndim, double *G, double *minusG);
int set_ynew(int ndim, double *yold, double *dely, double *ynew);
int set_yold(int ndim, double *yold, double *ynew);
int test_tol(int ndim, double *dely, double tol);

#ifdef TEST_NEWT_RAPH

int f1(int ndim, double *y, double *G, void *fdata);
int dfdy1(int ndim, double *y, double **dGdy, void *Jdata);
int f2(int ndim, double *y, double *G, void *fdata);
int dfdy2(int ndim, double *y, double **dGdy, void *Jdata);
int f3(int ndim, double *y, double *G, void *fdata);
int dfdy3(int ndim, double *y, double **dGdy, void *Jdata);

#endif

#define NR_HEADER_ALREADY_INCLUDED
#endif
