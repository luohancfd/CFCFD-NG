/*************************************************************************
 *                                                                       *
 * File: cvdemk.c                                                        *
 * Programmers: Scott D. Cohen and Alan C. Hindmarsh @ LLNL              *
 * Version of 25 February 2000                                           * 
 *-----------------------------------------------------------------------*
 * Demonstration program for CVODE - Krylov linear solver.               *
 * ODE system from ns-species interaction PDE in 2 dimensions.           *
 *                                                                       *
 * This program solves a stiff ODE system that arises from a system      *
 * of partial differential equations.  The PDE system is a food web      *
 * population model, with predator-prey interaction and diffusion on     *
 * the unit square in two dimensions.  The dependent variable vector is  *
 *                                                                       *
 *        1   2        ns                                                *
 *  c = (c , c , ..., c  )                                               *
 *                                                                       *
 * and the PDEs are as follows:                                          *
 *                                                                       *
 *    i               i      i                                           *
 *  dc /dt  =  d(i)*(c    + c   )  +  f (x,y,c)  (i=1,...,ns)            *
 *                    xx     yy        i                                 *
 *                                                                       *
 * where                                                                 *
 *                                                                       * 
 *                 i          ns         j                               *
 *  f (x,y,c)  =  c *(b(i) + sum a(i,j)*c )   .                          *
 *   i                       j=1                                         *
 *                                                                       *
 * The number of species is ns = 2*np, with the first np being prey and  *
 * the last np being predators.  The coefficients a(i,j), b(i), d(i) are *
 *                                                                       *
 *  a(i,i) = -a  (all i)                                                 *
 *  a(i,j) = -g  (i <= np, j > np)                                       *
 *  a(i,j) =  e  (i > np, j <= np)                                       *
 *  b(i) =  b*(1 + alpha*x*y)  (i <= np)                                 *
 *  b(i) = -b*(1 + alpha*x*y)  (i > np)                                  *
 *  d(i) = Dprey  (i <= np)                                              *
 *  d(i) = Dpred  (i > np)                                               *
 *                                                                       *
 * The spatial domain is the unit square. The final time is 10.          *
 * The boundary conditions are: normal derivative = 0.                   *
 * A polynomial in x and y is used to set the initial conditions.        *
 *                                                                       *
 * The PDEs are discretized by central differencing on an MX by MY mesh. *
 * The resulting ODE system is stiff.                                    *
 *                                                                       *
 * The ODE system is solved by CVODE using Newton iteration and the      *
 * CVSPGMR linear solver (scaled preconditioned GMRES).                  *
 *                                                                       *
 * The preconditioner matrix used is the product of two matrices:        * 
 * (1) A matrix, only defined implicitly, based on a fixed number of     *
 * Gauss-Seidel iterations using the diffusion terms only.               *
 * (2) A block-diagonal matrix based on the partial derivatives of the   *
 * interaction terms f only, using block-grouping (computing only a      * 
 * subset of the ns by ns blocks).                                       *
 *                                                                       *
 * Four different runs are made for this problem.                        *
 * The product preconditoner is applied on the left and on the right.    *
 * In each case, both the modified and classical Gram-Schmidt options    *
 * are tested.                                                           *
 * In the series of runs, CVodeMalloc is called only for the first run,  *
 * whereas CVReInit is called for each of the remaining seven runs.      *
 *                                                                       *
 * A problem description, performance statistics, and solution profiles  *
 * at selected output times are written to standard output. CVODE error  *
 * and warning messages are written to standard error. There should be   *
 * no such messages.                                                     *
 *                                                                       *
 * Note.. This program requires the "small" dense linear solver routines *
 * denalloc, denallocpiv, denaddI, gefa, gesl, denfreepiv and denfree.   *
 *                                                                       *
 * Note.. This program assumes the sequential implementation for the     *
 * type N_Vector and uses the N_VDATA macro to gain access to the        *
 * contiguous array of components of an N_Vector.                        *
 *-----------------------------------------------------------------------*
 *                                                                       *
 * Reference..  Peter N. Brown and Alan C. Hindmarsh, Reduced Storage    *
 * Matrix Methods in Stiff ODE Systems, J. Appl. Math. & Comp., 31       *
 * (1989), pp. 40-91.  Also available as Lawrence Livermore National     *
 * Laboratory Report UCRL-95088, Rev. 1, June 1987.                      *
 *************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "llnltyps.h"  /* contains the definition for real, integer, boole  */
#include "cvode.h"     /* main CVODE header file                            */
#include "iterativ.h"  /* for types of preconditioning and Gram-Schmidt     */
#include "cvspgmr.h"   /* use CVSPGMR linear solver each internal step      */
#include "dense.h"     /* use generic DENSE linear solver for "small" dense */
                       /* matrix blocks in right preconditioner             */
#include "nvector.h"   /* contains the definition of type N_Vector          */
#include "llnlmath.h"  /* contains UnitRoundoff, RSqrt, SQR functions       */

/* Problem Specification Constants */

#define AA 1.0       /* AA = a */
#define EE 1e4       /* EE = e */
#define GG 0.5e-6    /* GG = g */
#define BB 1.0       /* BB = b */
#define DPREY 1.0    
#define DPRED 0.5
#define ALPH 1.0
#define NP 3
#define NS (2*NP)

/* Method Constants */

#define MX   6
#define MY   6
#define MXNS (MX*NS)
#define AX   1.0
#define AY   1.0
#define DX   (AX/(real)(MX-1))
#define DY   (AY/(real)(MY-1))
#define MP   NS
#define MQ   (MX*MY)
#define MXMP  (MX*MP)
#define NGX   2
#define NGY   2
#define NGRP  (NGX*NGY)
#define ITMAX 5

/* CVodeMalloc Constants */

#define NEQ   (NS*MX*MY)
#define T0    0.0
#define RTOL  1e-5
#define ATOL  1e-5
#define LMM   BDF
#define ITER  NEWTON
#define ITOL  SS
#define ERRFP stderr
#define OPTIN FALSE

/* CVSpgmr Constants */

#define MAXL 0            /* => use default = MIN(NEQ, 5)            */
#define DELT 0.0          /* => use default = 0.05                   */

/* Output Constants */

#define T1        1e-8
#define TOUT_MULT 10.0
#define DTOUT     1.0
#define NOUT      18

/* Note: The value for species i at mesh point (j,k) is stored in */
/* component number (i-1) + j*NS + k*NS*MX of an N_Vector,        */
/* where 1 <= i <= NS, 0 <= j < MX, 0 <= k < MY.                  */

/* Structure for user data */

typedef struct {
  real   **P[NGRP];
  integer *pivot[NGRP];
  int ns,  mxns, mp, mq, mx, my, ngrp, ngx, ngy, mxmp;
  int jgx[NGX+1], jgy[NGY+1], jigx[MX], jigy[MY];
  int jxr[NGX], jyr[NGY];
  real acoef[NS][NS], bcoef[NS], diff[NS];
  real cox[NS], coy[NS], dx, dy, srur;
  real fsave[NEQ];
} *WebData;

/* Private Helper Functions */

static WebData AllocUserData(void);
static void InitUserData(WebData wdata);
static void SetGroups(int m, int ng, int jg[], int jig[], int jr[]);
static void CInit(N_Vector c, WebData wdata);
static void PrintIntro(void);
static void PrintAllSpecies(N_Vector c, int ns, int mxns, real t);
static void PrintOutput(real t, long int iopt[], real ropt[]);
static void PrintFinalStats(long int iopt[]);
static void FreeUserData(WebData wdata);
static void WebRates(real x, real y, real t, real c[], real rate[],
		     WebData wdata);
static void fblock (real t, real cdata[], int jx, int jy, real cdotdata[],
		    WebData wdata);
static void GSIter(int N, real gamma, N_Vector z, N_Vector x, WebData wdata);

/* Small Vector Kernels */

static void v_inc_by_prod(real u[], real v[], real w[], int n);
static void v_sum_prods(real u[], real p[], real q[], real v[], real w[],
			int n);
static void v_prod(real u[], real v[], real w[], int n);
static void v_zero(real u[], int n);

/* Functions Called By The CVODE Solver */

static void f(integer N, real t, N_Vector y, N_Vector ydot, void *f_data);

static int Precond(integer N, real tn, N_Vector c, N_Vector fc,
		   boole jok, boole *jcurPtr, real gamma, N_Vector ewt, real h,
		   real uround, long int *nfePtr, void *P_data,
		   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3);

static int PSolve(integer N, real tn, N_Vector c, N_Vector fc,
		  N_Vector vtemp, real gamma, N_Vector ewt, real delta,
		  long int *nfePtr, N_Vector r, int lr, void *P_data,
		  N_Vector z);
     
/* Implementation */

main()
{
  real ropt[OPT_SIZE], abstol=ATOL, reltol=RTOL, t, tout;
  long int iopt[OPT_SIZE]; 
  N_Vector c;
  WebData wdata;
  void *cvode_mem;
  boole firstrun;
  int jpre, gstype, iout, flag, ns, mxns;

  /* Initializations */
  c = N_VNew(NEQ, NULL);
  wdata = AllocUserData();
  InitUserData(wdata);

  /* Print problem description */
  PrintIntro();

  /* Loop over jpre and gstype (four cases) */
  for (jpre = LEFT; jpre <= RIGHT; jpre++) {
    for (gstype = MODIFIED_GS; gstype <= CLASSICAL_GS; gstype++) {
      
      /* Initialize c and print heading */
      CInit(c, wdata);
      printf("\n\nPreconditioner type is           jpre = %s\n",
	     (jpre == LEFT) ? "LEFT" : "RIGHT");
      printf("\n\nGram-Schmidt method type is    gstype = %s\n\n\n",
	     (gstype == MODIFIED_GS) ? "MODIFIED_GS" : "CLASSICAL_GS");

      /* Call CVodeMalloc or CVReInit, then CVSpgmr to set up problem */

      firstrun = (jpre==LEFT) && (gstype==MODIFIED_GS);
      if (firstrun) {
        cvode_mem = CVodeMalloc(NEQ, f, T0, c, LMM, ITER, ITOL, &reltol,
			  &abstol, wdata, ERRFP, OPTIN, iopt, ropt, NULL);
        if (cvode_mem == NULL) { printf("CVodeMalloc failed.\n"); return(1); }
      }
      else {
        flag = CVReInit(cvode_mem, f, T0, c, LMM, ITER, ITOL, &reltol,
			  &abstol, wdata, ERRFP, OPTIN, iopt, ropt, NULL);
        if (flag != SUCCESS) { printf("CVReInit failed."); return(1); }
      }

      CVSpgmr(cvode_mem, jpre, gstype, MAXL, DELT, Precond, PSolve, wdata);
 
      /* Print initial values */
      if ((jpre == LEFT) && (gstype == MODIFIED_GS))
         PrintAllSpecies(c, ns = wdata->ns, mxns = wdata->mxns, T0);

      /* Loop over output points, call CVode, print sample solution values. */
      tout = T1;
      for (iout = 1; iout <= NOUT; iout++) {
        flag = CVode(cvode_mem, tout, c, &t, NORMAL);
        PrintOutput(t, iopt, ropt);
        if ((jpre == LEFT) && (gstype == MODIFIED_GS) && (iout % 3 == 0)) 
                         PrintAllSpecies(c, ns, mxns, t);
        if (flag != SUCCESS) { printf("CVode failed.\n"); break; }
        if (tout > 0.9) tout += DTOUT; else tout *= TOUT_MULT; 
      }

      /* Print final statistics, and loop for next case */
      PrintFinalStats(iopt);

    }
  }

  /* Free all memory */
  CVodeFree(cvode_mem);
  N_VFree(c);
  FreeUserData(wdata);
  return(0);
}

static WebData AllocUserData(void)
{
  int i, ngrp = NGRP, ns = NS;
  WebData wdata;

  wdata = (WebData) malloc(sizeof *wdata);
  for(i=0; i < ngrp; i++) {
    (wdata->P)[i] = denalloc(ns);
    (wdata->pivot)[i] = denallocpiv(ns);
  }
  return(wdata);
}

static void InitUserData(WebData wdata)
{
  int i, j, ns;
  real *bcoef, *diff, *cox, *coy, dx, dy;
  real (*acoef)[NS];

  acoef = wdata->acoef;
  bcoef = wdata->bcoef;
  diff = wdata->diff;
  cox = wdata->cox;
  coy = wdata->coy;
  ns = wdata->ns = NS;

  for (j = 0; j < NP; j++) {
    for (i = 0; i < NP; i++) {
      acoef[NP+i][j] = EE;
      acoef[i][NP+j] = -GG;
    }
    acoef[j][j] = -AA;
    acoef[NP+j][NP+j] = -AA;
    bcoef[j] = BB;
    bcoef[NP+j] = -BB;
    diff[j] = DPREY;
    diff[NP+j] = DPRED;
  }

  /* Set remaining problem parameters */

  wdata->mxns = MXNS;
  dx = wdata->dx = DX;
  dy = wdata->dy = DY;
  for (i = 0; i < ns; i++) {
    cox[i] = diff[i]/SQR(dx);
    coy[i] = diff[i]/SQR(dy);
  }

  /* Set remaining method parameters */

  wdata->mp = MP;
  wdata->mq = MQ;
  wdata->mx = MX;
  wdata->my = MY;
  wdata->srur = RSqrt(UnitRoundoff());
  wdata->mxmp = MXMP;
  wdata->ngrp = NGRP;
  wdata->ngx = NGX;
  wdata->ngy = NGY;
  SetGroups(MX, NGX, wdata->jgx, wdata->jigx, wdata->jxr);
  SetGroups(MY, NGY, wdata->jgy, wdata->jigy, wdata->jyr);
}

/*
 This routine sets arrays jg, jig, and jr describing
 a uniform partition of (0,1,2,...,m-1) into ng groups.
 The arrays set are:
   jg    = length ng+1 array of group boundaries.
           Group ig has indices j = jg[ig],...,jg[ig+1]-1.
   jig   = length m array of group indices vs node index.
           Node index j is in group jig[j].
   jr    = length ng array of indices representing the groups.
           The index for group ig is j = jr[ig].
*/
static void SetGroups(int m, int ng, int jg[], int jig[], int jr[])
{
  int ig, j, len1, mper, ngm1;
  
  mper = m/ng; /* does integer division */
  for (ig=0; ig < ng; ig++) jg[ig] = ig*mper;
  jg[ng] = m;
  
  ngm1 = ng - 1;
  len1 = ngm1*mper;
  for (j = 0; j < len1; j++) jig[j] = j/mper;
  for (j = len1; j < m; j++) jig[j] = ngm1;

  for (ig = 0; ig < ngm1; ig++) jr[ig] = ((2*ig+1)*mper-1)/2;
  jr[ngm1] = (ngm1*mper+m-1)/2;
}

/* This routine computes and loads the vector of initial values. */
static void CInit(N_Vector c, WebData wdata)
{
  int i, ici, ioff, iyoff, jx, jy, ns, mxns;
  real argx, argy, x, y, dx, dy, x_factor, y_factor, *cdata;
  
  cdata = N_VDATA(c);
  ns = wdata->ns;
  mxns = wdata->mxns;
  dx = wdata->dx;
  dy = wdata->dy;

  x_factor = 4.0/SQR(AX);
  y_factor = 4.0/SQR(AY);
  for (jy = 0; jy < MY; jy++) {
    y = jy*dy;
    argy = SQR(y_factor*y*(AY-y)); 
    iyoff = mxns*jy;
    for (jx = 0; jx < MX; jx++) {
      x = jx*dx;
      argx = SQR(x_factor*x*(AX-x));
      ioff = iyoff + ns*jx;
      for (i = 1; i <= ns; i++) {
	ici = ioff + i-1;
	cdata[ici] = 10.0 + i*argx*argy;
      }
    }
  }
}

static void PrintIntro(void)
{
  printf("\n\nDemonstration program for CVODE - CVSPGMR linear solver\n\n");
  printf("Food web problem with ns species, ns = %d\n", NS);
  printf("Predator-prey interaction and diffusion on a 2-D square\n\n");
  printf("Matrix parameters.. a = %.2g   e = %.2g   g = %.2g\n",
	 AA, EE, GG);
  printf("b parameter = %.2g\n", BB);
  printf("Diffusion coefficients.. Dprey = %.2g   Dpred = %.2g\n",
	 DPREY, DPRED);
  printf("Rate parameter alpha = %.2g\n\n", ALPH);
  printf("Mesh dimensions (mx,my) are %d, %d.  ", MX, MY);
  printf("Total system size is neq = %d \n\n", NEQ);
  printf("Tolerances: itol = %s,  reltol = %.2g, abstol = %.2g \n\n",
         (ITOL == SS) ? "SS" : "SV", RTOL, ATOL);

  printf("Preconditioning uses a product of:\n");
  printf("  (1) Gauss-Seidel iterations with ");
  printf("itmax = %d iterations, and\n", ITMAX);
  printf("  (2) interaction-only block-diagonal matrix ");
  printf("with block-grouping\n");
  printf("  Number of diagonal block groups = ngrp = %d", NGRP);
  printf("  (ngx by ngy, ngx = %d, ngy = %d)\n", NGX, NGY);
  printf("\n\n--------------------------------------------------------------");
  printf("--------------\n");
}

static void PrintAllSpecies(N_Vector c, int ns, int mxns, real t)
{
  int i, jx, jy;
  real *cdata;

  cdata = N_VDATA(c);
  printf("c values at t = %g:\n\n", t);
  for (i=1; i <= ns; i++) {
    printf("Species %d\n", i);
    for (jy=MY-1; jy >= 0; jy--) {
      for (jx=0; jx < MX; jx++) {
	printf("%-10.6g", cdata[(i-1) + jx*ns + jy*mxns]);
      }
      printf("\n");
    }
    printf("\n");
  }
}

static void PrintOutput(real t, long int iopt[], real ropt[])
{
  printf("t = %10.2e  nst = %ld  nfe = %ld  nni = %ld", t, iopt[NST],
	 iopt[NFE], iopt[NNI]);
  printf("  qu = %ld  hu = %11.2e\n\n", iopt[QU], ropt[HU]);
}

static void PrintFinalStats(long int iopt[])
{
  long int nni, nli;
  real avdim;

  printf("\n\n Final statistics for this run:\n\n");
  printf(" CVode real workspace length           = %4ld \n", iopt[LENRW]);
  printf(" CVode integer workspace length        = %4ld \n", iopt[LENIW]);
  printf(" CVSPGMR real workspace length         = %4ld \n", iopt[SPGMR_LRW]);
  printf(" CVSPGMR integer workspace length      = %4ld \n", iopt[SPGMR_LIW]);
  printf(" Number of steps                       = %4ld \n", iopt[NST]);
  printf(" Number of f-s                         = %4ld \n", iopt[NFE]);
  printf(" Number of setups                      = %4ld \n", iopt[NSETUPS]);
  printf(" Number of nonlinear iterations        = %4ld \n", nni = iopt[NNI]);
  printf(" Number of linear iterations           = %4ld \n", 
	 nli = iopt[SPGMR_NLI]);
  printf(" Number of preconditioner evaluations  = %4ld \n", iopt[SPGMR_NPE]);
  printf(" Number of preconditioner solves       = %4ld \n", iopt[SPGMR_NPS]);
  printf(" Number of error test failures         = %4ld \n", iopt[NETF]);
  printf(" Number of nonlinear conv. failures    = %4ld \n", iopt[NCFN]);
  printf(" Number of linear convergence failures = %4ld \n", iopt[SPGMR_NCFL]);
  avdim = (nni > 0) ? ((real)nli)/((real)nni) : 0.0;
  printf(" Average Krylov subspace dimension     = %.3f \n", avdim);
  printf("\n\n--------------------------------------------------------------");
  printf("--------------\n");
  printf(    "--------------------------------------------------------------");
  printf("--------------\n");
}

static void FreeUserData(WebData wdata)
{
  int i, ngrp;

  ngrp = wdata->ngrp;
  for(i=0; i < ngrp; i++) {
    denfree((wdata->P)[i]);
    denfreepiv((wdata->pivot)[i]);
  }
  free(wdata);
}

/*
 This routine computes the right-hand side of the ODE system and
 returns it in cdot. The interaction rates are computed by calls to WebRates,
 and these are saved in fsave for use in preconditioning.
*/
static void f(integer N, real t, N_Vector c, N_Vector cdot, void *f_data)
{
  int i, ic, ici, idxl, idxu, idyl, idyu, iyoff, jx, jy, ns, mxns;
  real dcxli, dcxui, dcyli, dcyui, x, y, *cox, *coy, *fsave, dx, dy;
  real *cdata, *cdotdata;
  WebData wdata;

  wdata = (WebData) f_data;
  cdata = N_VDATA(c);
  cdotdata = N_VDATA(cdot);

  mxns = wdata->mxns;
  ns = wdata->ns;
  fsave = wdata->fsave;
  cox = wdata->cox;
  coy = wdata->coy;
  mxns = wdata->mxns;
  dx = wdata->dx;
  dy = wdata->dy;

  for (jy = 0; jy < MY; jy++) {
    y = jy*dy;
    iyoff = mxns*jy;
    idyu = (jy == MY-1) ? -mxns : mxns;
    idyl = (jy == 0) ? -mxns : mxns;
    for (jx = 0; jx < MX; jx++) {
      x = jx*dx;
      ic = iyoff + ns*jx;
      /* Get interaction rates at one point (x,y). */
      WebRates(x, y, t, cdata+ic, fsave+ic, wdata);
      idxu = (jx == MX-1) ? -ns : ns;
      idxl = (jx == 0) ? -ns : ns;
      for (i = 1; i <= ns; i++) {
	ici = ic + i-1;
	/* Do differencing in y. */
	dcyli = cdata[ici] - cdata[ici-idyl];
	dcyui = cdata[ici+idyu] - cdata[ici];
	/* Do differencing in x. */
	dcxli = cdata[ici] - cdata[ici-idxl];
	dcxui = cdata[ici+idxu] - cdata[ici];
	/* Collect terms and load cdot elements. */
	cdotdata[ici] = coy[i-1]*(dcyui - dcyli) + cox[i-1]*(dcxui - dcxli) +
	  fsave[ici];
      }
    }
  }
}

/*
 This routine computes the interaction rates for the species
 c_1, ... ,c_ns (stored in c[0],...,c[ns-1]), at one spatial point 
 and at time t.
*/
static void WebRates(real x, real y, real t, real c[], real rate[],
		     WebData wdata)
{
  int i, j, ns;
  real fac, *bcoef;
  real (*acoef)[NS];

  ns = wdata->ns;
  acoef = wdata->acoef;
  bcoef = wdata->bcoef;

  for (i = 0; i < ns; i++)
    rate[i] = 0.0;
  
  for (j = 0; j < ns; j++) 
    for (i = 0; i < ns; i++) 
      rate[i] += c[j] * acoef[i][j];
  
  fac = 1.0 + ALPH*x*y;
  for (i = 0; i < ns; i++) 
    rate[i] = c[i]*(bcoef[i]*fac + rate[i]);
}

/*
 This routine generates the block-diagonal part of the Jacobian
 corresponding to the interaction rates, multiplies by -gamma, adds
 the identity matrix, and calls gefa to do the LU decomposition of
 each diagonal block. The computation of the diagonal blocks uses
 the preset block and grouping information. One block per group is
 computed. The Jacobian elements are generated by difference
 quotients using calls to the routine fblock.

 This routine can be regarded as a prototype for the general case
 of a block-diagonal preconditioner. The blocks are of size mp, and
 there are ngrp=ngx*ngy blocks computed in the block-grouping scheme.
*/ 
static int Precond(integer N, real t, N_Vector c, N_Vector fc, boole jok,
		   boole *jcurPtr, real gamma, N_Vector rewt, real h,
		   real uround, long int *nfePtr, void *P_data,
		   N_Vector vtemp1, N_Vector vtemp2, N_Vector vtemp3)
{
  real ***P;
  integer **pivot, ier;
  int i, if0, if00, ig, igx, igy, j, jj, jx, jy;
  int *jxr, *jyr, mp, ngrp, ngx, ngy, mxmp;
  real fac, r, r0, save, srur;
  real *f1, *fsave, *cdata, *rewtdata;
  WebData wdata;

  wdata = (WebData) P_data;
  cdata = N_VDATA(c);
  rewtdata = N_VDATA(rewt);

  P = wdata->P;
  pivot = wdata->pivot;
  jxr = wdata->jxr;
  jyr = wdata->jyr;
  mp = wdata->mp;
  srur = wdata->srur;
  ngrp = wdata->ngrp;
  ngx = wdata->ngx;
  ngy = wdata->ngy;
  mxmp = wdata->mxmp;
  fsave = wdata->fsave;

  /* Make mp calls to fblock to approximate each diagonal block of Jacobian.
     Here, fsave contains the base value of the rate vector and 
     r0 is a minimum increment factor for the difference quotient. */

  f1 = N_VDATA(vtemp1);

  fac = N_VWrmsNorm (fc, rewt);
  r0 = 1000.0*ABS(gamma)*uround*N*fac;
  if (r0 == 0.0) r0 = 1.0;

  for (igy = 0; igy < ngy; igy++) {
    jy = jyr[igy];
    if00 = jy*mxmp;
    for (igx = 0; igx < ngx; igx++) { 
      jx = jxr[igx];
      if0 = if00 + jx*mp;
      ig = igx + igy*ngx; 
      /* Generate ig-th diagonal block */
      for (j = 0; j < mp; j++) {
	/* Generate the jth column as a difference quotient */
	jj = if0 + j; 
	save = cdata[jj];
	r = MAX(srur*ABS(save),r0/rewtdata[jj]);
	cdata[jj] += r;
	fac = -gamma/r;
	fblock (t, cdata, jx, jy, f1, wdata);
	for (i = 0; i < mp; i++) {
	  P[ig][j][i] = (f1[i] - fsave[if0+i])*fac;
	}
	cdata[jj] = save;
      }
    }
  }

  /* Add identity matrix and do LU decompositions on blocks. */

   for (ig = 0; ig < ngrp; ig++) {
     denaddI(P[ig], mp);
     ier = gefa(P[ig], mp, pivot[ig]);
     if (ier != 0) return(1);
   }

  *jcurPtr = TRUE;
  return(0);
}

/*
 This routine computes one block of the interaction terms of the
 system, namely block (jx,jy), for use in preconditioning.
 Here jx and jy count from 0.
*/
static void fblock(real t, real cdata[], int jx, int jy, real cdotdata[],
		   WebData wdata)
{
  int iblok, ic;
  real x, y;

  iblok = jx + jy*(wdata->mx);
  y = jy*(wdata->dy);
  x = jx*(wdata->dx);
  ic = (wdata->ns)*(iblok);
  WebRates(x, y, t, cdata+ic, cdotdata, wdata);
}

/*
 This routine applies two inverse preconditioner matrices
 to the vector r, using the interaction-only block-diagonal Jacobian
 with block-grouping, denoted Jr, and Gauss-Seidel applied to the
 diffusion contribution to the Jacobian, denoted Jd.
 It first calls GSIter for a Gauss-Seidel approximation to
 ((I - gamma*Jd)-inverse)*r, and stores the result in z.
 Then it computes ((I - gamma*Jr)-inverse)*z, using LU factors of the
 blocks in P, and pivot information in pivot, and returns the result in z.
*/
static int PSolve(integer N, real tn, N_Vector c, N_Vector fc, N_Vector vtemp,
		  real gamma, N_Vector rewt, real delta, long int *nfePtr,
		  N_Vector r, int lr, void *P_data, N_Vector z)
{
  real   ***P;
  integer **pivot;
  int jx, jy, igx, igy, iv, ig, *jigx, *jigy, mx, my, ngx, mp;
  WebData wdata;

  wdata = (WebData) P_data;

  N_VScale(1.0, r, z);

  /* call GSIter for Gauss-Seidel iterations */

  GSIter(N, gamma, z, vtemp, wdata);

  /* Do backsolves for inverse of block-diagonal preconditioner factor */
 
  P = wdata->P;
  pivot = wdata->pivot;
  mx = wdata->mx;
  my = wdata->my;
  ngx = wdata->ngx;
  mp = wdata->mp;
  jigx = wdata->jigx;
  jigy = wdata->jigy;

  iv = 0;
  for (jy = 0; jy < my; jy++) {
    igy = jigy[jy];
    for (jx = 0; jx < mx; jx++) {
      igx = jigx[jx];
      ig = igx + igy*ngx;
      gesl(P[ig], mp, pivot[ig], &(N_VDATA(z)[iv]));
      iv += mp;
    }
  }

  return(0);
}

/*
 This routine performs ITMAX=5 Gauss-Seidel iterations to compute an
 approximation to (P-inverse)*z, where P = I - gamma*Jd, and
 Jd represents the diffusion contributions to the Jacobian.
 The answer is stored in z on return, and x is a temporary vector.
 The dimensions below assume a global constant NS >= ns.
 Some inner loops of length ns are implemented with the small
 vector kernels v_sum_prods, v_prod, v_inc_by_prod.
*/
static void GSIter(int N, real gamma, N_Vector z, N_Vector x, WebData wdata)
{
  int i, ic, ici, iter, iyoff, jx, jy, ns, mxns, mx, my, x_loc, y_loc;
  real beta[NS], beta2[NS], cof1[NS], gam[NS], gam2[NS];
  real temp, *cox, *coy, *xd, *zd;

  xd = N_VDATA(x);
  zd = N_VDATA(z);
  ns = wdata->ns;
  mx = wdata->mx;
  my = wdata->my;
  mxns = wdata->mxns;
  cox = wdata->cox;
  coy = wdata->coy;

  /* Write matrix as P = D - L - U.
     Load local arrays beta, beta2, gam, gam2, and cof1. */
 
  for (i = 0; i < ns; i++) {
    temp = 1.0/(1.0 + 2.0*gamma*(cox[i] + coy[i]));
    beta[i] = gamma*cox[i]*temp;
    beta2[i] = 2.0*beta[i];
    gam[i] = gamma*coy[i]*temp;
    gam2[i] = 2.0*gam[i];
    cof1[i] = temp;
  }

  /* Begin iteration loop.
  Load vector x with (D-inverse)*z for first iteration. */

  for (jy = 0; jy < my; jy++) {
    iyoff = mxns*jy;
    for (jx = 0; jx < mx; jx++) {
      ic = iyoff + ns*jx;
      v_prod(xd+ic, cof1, zd+ic, ns); /* x[ic+i] = cof1[i]z[ic+i] */
    }
  }
  N_VConst(0.0, z);

  /* Looping point for iterations. */

  for (iter=1; iter <= ITMAX; iter++) {

    /* Calculate (D-inverse)*U*x if not the first iteration. */

    if (iter > 1) {
      for (jy=0; jy < my; jy++) {
	iyoff = mxns*jy;
	for (jx=0; jx < mx; jx++) { /* order of loops matters */
	  ic = iyoff + ns*jx;
	  x_loc = (jx == 0) ? 0 : ((jx == mx-1) ? 2 : 1);
	  y_loc = (jy == 0) ? 0 : ((jy == my-1) ? 2 : 1);
	  switch (3*y_loc+x_loc) {
	  case 0 : /* jx == 0, jy == 0 */
	           /* x[ic+i] = beta2[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] */
	           v_sum_prods(xd+ic, beta2, xd+ic+ns, gam2, xd+ic+mxns, ns);
		   break;
	  case 1 : /* 1 <= jx <= mx-2, jy == 0 */
	           /* x[ic+i] = beta[i]x[ic+ns+i] + gam2[i]x[ic+mxns+i] */
	           v_sum_prods(xd+ic, beta, xd+ic+ns, gam2, xd+ic+mxns, ns);
		   break;
	  case 2 : /* jx == mx-1, jy == 0 */
	           /* x[ic+i] = gam2[i]x[ic+mxns+i] */
	           v_prod(xd+ic, gam2, xd+ic+mxns, ns);
		   break;
	  case 3 : /* jx == 0, 1 <= jy <= my-2 */
	           /* x[ic+i] = beta2[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] */
                   v_sum_prods(xd+ic, beta2, xd+ic+ns, gam, xd+ic+mxns, ns);
		   break;
	  case 4 : /* 1 <= jx <= mx-2, 1 <= jy <= my-2 */
	           /* x[ic+i] = beta[i]x[ic+ns+i] + gam[i]x[ic+mxns+i] */
	           v_sum_prods(xd+ic, beta, xd+ic+ns, gam, xd+ic+mxns, ns);
		   break;
	  case 5 : /* jx == mx-1, 1 <= jy <= my-2 */
	           /* x[ic+i] = gam[i]x[ic+mxns+i] */
	           v_prod(xd+ic, gam, xd+ic+mxns, ns);
		   break;
	  case 6 : /* jx == 0, jy == my-1 */
	           /* x[ic+i] = beta2[i]x[ic+ns+i] */
	           v_prod(xd+ic, beta2, xd+ic+ns, ns);
		   break;
	  case 7 : /* 1 <= jx <= mx-2, jy == my-1 */
	           /* x[ic+i] = beta[i]x[ic+ns+i] */
	           v_prod(xd+ic, beta, xd+ic+ns, ns);
		   break;
	  case 8 : /* jx == mx-1, jy == my-1 */
	           /* x[ic+i] = 0.0 */
	           v_zero(xd+ic, ns);
		   break;
	  }
	}
      }
    }  /* end if (iter > 1) */
 
    /* Overwrite x with [(I - (D-inverse)*L)-inverse]*x. */

    for (jy=0; jy < my; jy++) {
      iyoff = mxns*jy;
      for (jx=0; jx < mx; jx++) { /* order of loops matters */
	ic = iyoff + ns*jx;
	x_loc = (jx == 0) ? 0 : ((jx == mx-1) ? 2 : 1);
	y_loc = (jy == 0) ? 0 : ((jy == my-1) ? 2 : 1);
	switch (3*y_loc+x_loc) {
	case 0 : /* jx == 0, jy == 0 */
	         break;
	case 1 : /* 1 <= jx <= mx-2, jy == 0 */
	         /* x[ic+i] += beta[i]x[ic-ns+i] */
	         v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
		 break;
	case 2 : /* jx == mx-1, jy == 0 */
	         /* x[ic+i] += beta2[i]x[ic-ns+i] */
	         v_inc_by_prod(xd+ic, beta2, xd+ic-ns, ns);
		 break;
	case 3 : /* jx == 0, 1 <= jy <= my-2 */
	         /* x[ic+i] += gam[i]x[ic-mxns+i] */
	         v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
		 break;
	case 4 : /* 1 <= jx <= mx-2, 1 <= jy <= my-2 */
	         /* x[ic+i] += beta[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] */
	         v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
		 v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
		 break;
	case 5 : /* jx == mx-1, 1 <= jy <= my-2 */
	         /* x[ic+i] += beta2[i]x[ic-ns+i] + gam[i]x[ic-mxns+i] */
	         v_inc_by_prod(xd+ic, beta2, xd+ic-ns, ns);
		 v_inc_by_prod(xd+ic, gam, xd+ic-mxns, ns);
		 break;
	case 6 : /* jx == 0, jy == my-1 */
	         /* x[ic+i] += gam2[i]x[ic-mxns+i] */
	         v_inc_by_prod(xd+ic, gam2, xd+ic-mxns, ns);
		 break;
	case 7 : /* 1 <= jx <= mx-2, jy == my-1 */
	         /* x[ic+i] += beta[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] */
	         v_inc_by_prod(xd+ic, beta, xd+ic-ns, ns);
		 v_inc_by_prod(xd+ic, gam2, xd+ic-mxns, ns);
                 break;
	case 8 : /* jx == mx-1, jy == my-1 */
	         /* x[ic+i] += beta2[i]x[ic-ns+i] + gam2[i]x[ic-mxns+i] */
	         v_inc_by_prod(xd+ic, beta2, xd+ic-ns, ns);
		 v_inc_by_prod(xd+ic, gam2, xd+ic-mxns, ns);
		 break;
	}
      }
    }
 
    /* Add increment x to z : z <- z+x */
    
    N_VLinearSum(1.0, z, 1.0, x, z);

  }
}

static void v_inc_by_prod(real u[], real v[], real w[], int n)
{
  int i;  
  for (i=0; i < n; i++) u[i] += v[i]*w[i];
}

static void v_sum_prods(real u[], real p[], real q[], real v[], real w[],
			int n)
{
  int i;  
  for (i=0; i < n; i++) u[i] = p[i]*q[i] + v[i]*w[i];
}

static void v_prod(real u[], real v[], real w[], int n)
{ 
  int i;
  for (i=0; i < n; i++) u[i] = v[i]*w[i];
}

static void v_zero(real u[], int n)
{
  int i;  
  for (i=0; i < n; i++) u[i] = 0.0;
}
