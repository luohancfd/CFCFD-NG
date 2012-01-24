/** \file vecmat.h
 * \ingroup nm
 * \brief Header file for vecmat.c
 *
 * \date 24-Apr-04
 * \author Rowan J Gollan
 * \version 1.0
 *
 */
#ifndef VM_HEADER_ALREADY_INCLUDED

double *mem_alloc_vector(int ndim);
int *mem_alloc_iarray(int ndim);
double **mem_alloc_sqmat(int ndim);
int free_sqmat(double **A, int ndim);
int copy_vector(double *a, double *b, int ndim);
int print_vector(int ndim, double *vec, char *name);
int scale_vector(int ndim, double factor, double *v_in, double *v_out);
int add_vectors(int ndim, double *a, double *b, double *c);
int subtract_vectors(int ndim, double *a, double *b, double *c);
int add_scalar_to_vector(int ndim, double scalar, double *a, double *b);
int element_divide(int ndim, double *a, double *b, double *c);
int abs_vector(int ndim, double *a, double *b);
double max_vector(int ndim, double *a);
int print_sqmat(int ndim, double **sqmat, char *name);
int eye(double **I, int ndim);
int msys_subtract(double **A, double **B, double **C, int ndim);
int msys_scalmult(double C, double **A, double **B, int ndim);
int msys_solve(double **A, double *x, double *b, int ndim);


#define VM_HEADER_ALREADY_INCLUDED
#endif
